# !/usr/local/bin/python3
#
# transportation_problem.py
#
# 用单纯形算法求解运输问题
# 参考 [PulP Docs -- A Transportation Problem](https://coin-or.github.io/pulp/CaseStudies/a_transportation_problem.html)

import numpy as np
from linprog.problem import LpProblem
from linprog.simplex_vectorized import SimplexVectorized
from linprog.simplex_method import SimplexMethod


class TransportationProblem(object):
    class TransportationResult(object):
        def __init__(self, transportation_dict, total_cost):
            self.transportation_dict = transportation_dict
            self.total_cost = total_cost

    def __init__(self, supply: list, demand: list, costs: list):
        super().__init__()
        self.supply = supply
        self.demand = demand
        self.costs = costs
        self.result = None

    def solve(self, print_result=True) -> TransportationResult:
        # 将运输问题转化为线性规划问题
        lp_c = [-i for v in self.costs for i in v]  # 注意 `-i`，因为 linprog.problem 中 LpProblem 表示的是 max，而这里运输问题是 min。
        lp_a = []
        lp_b = []
        rows, cols = len(self.costs), len(self.costs[0])
        for i in range(len(self.supply)):
            a = [0] * len(lp_c)
            a[cols * i: cols * (i + 1)] = [1] * cols
            lp_a.append(a)
            lp_b.append(self.supply[i][1])
        for i in range(len(self.demand)):
            a = [1 if j % cols == i else 0 for j in range(len(lp_c))]
            lp_a.append(a)
            lp_b.append(self.demand[i][1])

        lpPb = LpProblem(lp_c, lp_a, lp_b)

        # 用单纯形法求解线性规划问题
        s = lpPb.solve(SimplexVectorized)
        if s.success:
            self.result = TransportationProblem.TransportationResult(np.reshape(s.solve, (rows, cols)), s.target * -1)

        # 输出解
        if print_result:
            if s.success:
                transportation = np.reshape(s.solve, (rows, cols))
                echo = [[' '] + [i[0] for i in self.demand]]
                for i in range(rows):
                    echo.append([self.supply[i][0]] + [j for j in transportation[i]])
                es = ''
                for r in echo:
                    for c in r:
                        es += f'{c}\t'
                    es += '\n'
                print(es)
            else:
                print("求解失败")

        return self.result


if __name__ == "__main__":
    supply = [('A1', 14), ('A2', 27), ('A3', 19)]
    demand = [('B1', 22), ('B2', 13), ('B3', 12), ('B4', 13)]
    costs = [[6, 7, 5, 3], [8, 4, 2, 7], [5, 9, 10, 6]]
    p = TransportationProblem(supply, demand, costs)
    p.solve()
