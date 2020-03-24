from linprog.problem import LpProblem
from linprog.simplex_method import SimplexMethod

if __name__ == "__main__":
    # 普通标准问题
    # pb = LpProblem([-5, 5, 13, 0, 0], [[-1, 1, 3, 1, 0], [12, 4, 10, 0, 1]], [20, 90])
    # s = pb.solve(SimplexMethod, show_tab=True)
    # print(s)

    # 唯一解
    pb = LpProblem([2, 3, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    s = pb.solve(SimplexMethod, show_tab=True)
    print(s)

    # 无穷多解
    # pb = LpProblem([2, 4, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    # s = pb.solve(SimplexMethod, show_tab=True)
    # print(s)

    # 无界解
    # pb = LpProblem([1, 1, 0, 0], [[-2, 1, 1, 0], [1, -1, 0, 1]], [4, 2])
    # s = pb.solve(SimplexMethod, show_tab=True)
    # print(s)

    # 大 M 法
    # pb = LpProblem([-2, -3, -1, 0, 0], [[1, 4, 2, -1, 0], [3, 2, 0, 0, -1]], [8, 6])
    # s = pb.solve(SimplexMethod, show_tab=True, big_m=True)
    # print(s)

    # 两阶段法
    # pb = LpProblem([-2, -3, -1, 0, 0], [[1, 4, 2, -1, 0], [3, 2, 0, 0, -1]], [8, 6])
    # s = pb.solve(SimplexMethod, show_tab=True, two_step=True)
    # print(s)

# import copy
#
# import numpy as np
#
#
# class LpProblem(object):
#     """
#     LpProblem 描述一个线性规划问题
#
#     :attribute c: 目标函数系数,  <math> c = (c_1, c_2, ..., c_n)^T </math>
#     :attribute a: 系数矩阵,     <math> a = (a_{ij})_{m \times n} = (p_1, p_2, ..., p_n) </math>
#     :attribute b: 右端常数,     <math> b = (b_1, b_2, ..., b_m)^T, b > 0 </math>
#     :attribute base_idx: 基变量的下标集合
#     """
#
#     def __init__(self, c: np.array, a: np.array, b: np.array):
#         self.c = np.array(c, 'float64')
#         self.a = np.array(a, 'float64')
#         self.b = np.array(b, 'float64')
#
#         self.base_idx = []
#
#         # For SimplexMethod
#         self.entering_idx = -1
#         self.leaving_idx = -1
#         self.theta = []
#         self.tab = []
#
#
# class LpSolve(object):
#     """
#     LpSolve 描述一个线性规划问题的解
#
#     :attributes success:     是否得到了最优解
#     :attributes description: 解的描述
#     :attributes solve:       最优解
#     """
#
#     def __init__(self, success: bool, description: str, solve: np.array):
#         self.success = success
#         self.description = description
#         self.solve = solve
#
#
# class SimplexMethod(object):
#     pass
#
#
# def simplex_solve(problem: LpProblem, tab=False) -> LpSolve:
#     """ simplex_solve 对给定<em>标准型</em>线性规划问题使用<em>单纯形</em>算法进行求解
#
#     待解决应符合标准型，即：
#         <math>
#             max z = c^T * x
#             s.t. a*x = b, x >= 0
#         </math>
#
#     单纯形算法参考：https://zh.wikipedia.org/zh-hans/单纯形法
#
#     可选参数：
#         tab=True (default False): 计算单纯形表
#
#     :return: 问题的解
#     """
#
#     p = copy.deepcopy(problem)
#
#     _current_simplex_tab(p)
#
#     p.entering_idx = np.argmax(p.c)  # 入基变量
#
#     while p.c[p.entering_idx] > 0:
#         p.theta = []
#         for i in range(len(p.b)):
#             if p.a[i][p.entering_idx] > 0:
#                 p.theta.append(p.b[i] / p.a[i][p.entering_idx])
#             else:
#                 p.theta.append(float("inf"))
#
#         p.leaving_idx = np.argmin(np.array(p.theta))  # 出基变量
#
#         if p.theta[p.leaving_idx] == float("inf"):  # 出基变量 == inf
#             return LpSolve(False, "出基变量 == inf", None)
#
#         _pivot(p)
#
#         _current_simplex_tab(p)
#
#         p.entering_idx = np.argmax(p.c)
#
#     x = np.zeros(len(p.c))
#     x[p.base_idx] = p.b
#     print(simplex_tab_tostring(p, None))
#     return LpSolve(True, "得到最优解", x)
#
#
# def _pivot(p: LpProblem):
#     """
#     对给定问题原址执行<em>转轴操作</em>（基变换）
#     """
#
#     main_element = p.a[p.leaving_idx][p.entering_idx]
#
#     p.a[p.leaving_idx] /= main_element
#     p.b[p.leaving_idx] /= main_element
#
#     p.base_idx[p.leaving_idx] = p.entering_idx
#
#     for i in range(len(p.b)):
#         if i != p.leaving_idx and p.a[i][p.entering_idx] != 0:
#             p.b[i] -= p.a[i][p.entering_idx] * p.b[p.leaving_idx]
#             p.a[i] -= p.a[i][p.entering_idx] * p.a[p.leaving_idx]
#
#     p.c -= p.c[p.entering_idx] * p.a[p.leaving_idx]
#
#
# def _current_simplex_tab(p: LpProblem):
#     """
#     计算当前单纯形表
#     :return: None
#     """
#     tab = [
#         (" ", " ") + tuple(p.c) + ("theta",),
#     ]
#     if len(p.tab) > 0:
#         main_element = '%.2f' % p.tab[-1][p.leaving_idx+1][p.entering_idx+3]
#         p.tab[-1][p.leaving_idx+1][p.entering_idx+2] = f'[{main_element}]'
#
#     for i in range(len(p.b)):
#         if len(p.tab) > 0 and len(p.theta) > 0:
#             p.tab[-1][i+1][-1] = p.theta[i]
#
#         tab += [
#             [f'x_{p.base_idx[i]}', p.b[i]] + list(p.a[i]) + [" ", ]
#         ]
#
#     p.tab.append(tab)
#
#
# def simplex_tab_tostring(p: LpProblem, step=None):
#     s = ''
#     if step is None:
#         for step in p.tab:
#             for row in step:
#                 for i in row:
#                     s += '%6.6s\t' % i
#                 s += '\n'
#             s += '-' * 16 + '\n'
#     else:
#         for row in p.tab[step]:
#             for i in row:
#                 s += '%6.6s\t' % i
#                 s += '\n'
#             s += '-' * 16 + '\n'
#     return s
#
#
# if __name__ == "__main__":
#     pb = LpProblem([-5, 5, 13, 0, 0], [[-1, 1, 3, 1, 0], [12, 4, 10, 0, 1]], [20, 90])
#     pb.base_idx = [3, 4]
#     s = simplex_solve(pb)
#     print(s.success)
#     print(s.description)
#     print(s.solve)
#     """
#     API Design:
#     pb = LpProblem([-5, 5, 13, 0, 0], [[-1, 1, 3, 1, 0], [12, 4, 10, 0, 1]], [20, 90])
#     s = pb.solve(Method, **kwargs)
#     """
