import copy

import numpy as np

from .problem import LpProblem
from .solve import LpSolve
from .solver import LpSolver


# TODO: Refactor: Abstract From SimplexVectorized and SimplexMethod, 面向切片完成单纯形表。
# TODO: Implement a faster inv() for base_inv


class SimplexVectorized(LpSolver):
    """ 矩阵化的单纯形法

    待解决应符合标准型，即：
    <math>
        max z = c^T * x
        s.t. a*x = b, x >= 0, b > 0
    </math>

    单纯形算法参考：https://zh.wikipedia.org/zh-hans/单纯形法
    """

    class Problem(LpProblem):
        """
        矩阵化的单纯形法内部的线性规划问题表示
        """
        def __init__(self, c, a, b):
            super().__init__(c, a, b)
            self.base_idx = np.ones(len(b), 'int') * -1
            self.cc = copy.deepcopy(c)

    def __init__(self, problem: LpProblem):
        super().__init__(problem)
        self.problem = self.Problem(problem.c, problem.a, problem.b)

    def find_idt_base(self):
        """
        尝试找「单位阵」作初始基

        :return: True if success False else
        """
        base_idx = np.ones(len(self.problem.b), 'int') * -1
        aT = self.problem.a.T
        for i in range(len(self.problem.b)):
            e = np.zeros(len(self.problem.b))
            e[i] = 1
            for j in range(len(aT)):
                if np.all(aT[j] == e):
                    base_idx[i] = j

        self.problem.base_idx = base_idx
        return np.all(base_idx >= 0)

    def big_m(self, **kwargs):
        """
        用大M法得到初始基

        :return: None
        """
        M = ((max(abs(self.problem.c)) + 1) ** 2) * 10 + 10

        for i in range(len(self.problem.base_idx)):
            if self.problem.base_idx[i] < 0:
                self.problem.c = np.insert(self.problem.c, len(self.problem.c), np.array([-M]))

                ap = np.zeros(len(self.problem.b))
                ap[i] = 1
                self.problem.a = np.c_[self.problem.a, ap]

                self.problem.base_idx[i] = len(self.problem.c) - 1

    def two_step(self, **kwargs):
        """
        用两阶段法得到初始基，第一阶段在此计算

        :return: 第一阶段的解
        """
        p = copy.deepcopy(self.problem)
        p.c = np.zeros(len(p.c))
        for i in range(len(self.problem.base_idx)):
            if self.problem.base_idx[i] < 0:
                p.c = np.insert(p.c, len(p.c), np.array([-1]))

                ap = np.zeros(len(p.b))
                ap[i] = 1
                p.a = np.c_[p.a, ap]

                self.problem.base_idx[i] = len(p.c) - 1
        p.base_idx = self.problem.base_idx
        s1 = _simplex_vectorized_solve(p)

        self.problem.c = copy.deepcopy(self.problem.cc)
        self.problem.base_idx = p.base_idx
        self.problem.b = p.b
        self.problem.a = (p.a.T[0: len(self.problem.c)]).T
        return s1

    def solve(self, **kwargs) -> LpSolve:
        """
        单纯形算法入口

        :param kwargs:
            base_idx=[...] (default []): 指定初始基，缺省则由算法自行确定
            two_step=True  (default False): 使用两阶段法
            big_m=True     (default True):  使用大 M 法
        :return: 问题的解
        """
        base_idx = kwargs.get("base_idx", [])
        if base_idx:    # 用户指定了基
            self.problem.base_idx = base_idx
        else:
            if not self.find_idt_base():    # 没有找到单位阵作初始基，用人工变量法（大M法 / 两阶段法）
                if kwargs.get("two_step", False):
                    s1 = self.two_step(**kwargs)
                    if not s1.success:  # 第一阶段确定无可行解
                        return s1
                else:
                    self.big_m(**kwargs)
        self.problem.base = np.identity(len(self.problem.b))
        self.problem.base_inv = np.identity(len(self.problem.b))

        pb = copy.deepcopy(self.problem)
        return _simplex_vectorized_solve(pb)


def _simplex_vectorized_solve(p: SimplexVectorized.Problem) -> LpSolve:
    """ _simplex_vectorized_solve 对给定了初始基的标准型使用<em>矩阵表示的单纯形法</em>进行求解

    :return: 问题的解
    """
    base_idx = p.base_idx
    n_base_idx = _get_n_base_idx(p, base_idx)

    base_inv = np.linalg.inv(p.a.T[base_idx].T)

    # 初始的检验数计算
    sigma_n = p.c[n_base_idx] - np.dot(np.dot(p.c[base_idx], base_inv), p.a.T[n_base_idx].T)

    while np.any(sigma_n > 0):
        # 确定入基变量
        e = np.argmax(sigma_n)
        entering_idx = n_base_idx[e]

        base_inv_p_entering = np.dot(base_inv, p.a.T[entering_idx].T)

        if np.all(base_inv_p_entering <= 0):
            return LpSolve(False, "无界解", None, None)

        # 确定换出变量
        theta = _get_theta(base_inv, base_inv_p_entering, p)

        l = np.argmin(theta)
        leaving_idx = base_idx[l]

        # 基变换
        base_idx[l] = entering_idx
        n_base_idx[e] = leaving_idx

        base_inv = np.linalg.inv(p.a.T[base_idx].T)
        # Next Sigma_N
        sigma_n = p.c[n_base_idx] - np.dot(np.dot(p.c[base_idx], base_inv), p.a.T[n_base_idx].T)

    # 迭代结束，计算解的值
    # base_inv = np.linalg.inv(p.a.T[base_idx].T)
    xb = np.dot(base_inv, p.b)
    x = np.zeros(len(p.c))
    x[base_idx] = xb

    # 分析解的情况
    return _analyse_solve(p, x, sigma_n)


def _get_theta(base_inv, base_inv_p_entering, p):
    """
    检验数计算

    :return: 检验数
    """
    theta = np.ones(len(p.b)) * np.inf
    for i in range(len(theta)):
        if base_inv_p_entering[i] > 0:
            theta[i] = np.dot(base_inv, p.b)[i] / base_inv_p_entering[i]
    return theta


def _get_n_base_idx(p: SimplexVectorized.Problem, base_idx):
    """
    初始化非基

    :param p: SimplexVectorized.Problem
    :param base_idx: 已确定的基indices
    :return: 非基indices
    """
    n_base_idx = list(range(len(p.c)))
    for i in base_idx:
        try:
            n_base_idx.pop(n_base_idx.index(i))
        except ValueError:
            pass
    return np.array(n_base_idx, 'int')


def _analyse_solve(p: SimplexVectorized.Problem, x, last_sigma_n) -> LpSolve:
    """
    单纯形求解完成后，分析解的情况

    :return: 单纯形求解结果(LpSolve实例)
    """
    x_real = x[0: len(p.cc)]
    x_personal = x[len(p.cc):]    # 人工变量

    if np.any(x_personal != 0):
        return LpSolve(False, "无可行解", None, None)

    z = np.dot(x_real, p.cc)

    for i in last_sigma_n:
        if abs(i) < 1e-8:  # 非基变量检验数为0
            return LpSolve(True, "无穷多最优解", x_real, z)

    return LpSolve(True, "唯一最优解", x_real, z)
