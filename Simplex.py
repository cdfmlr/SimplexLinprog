from linprog.problem import LpProblem
from linprog.simplex_method import SimplexMethod
from linprog.simplex_vectorized import SimplexVectorized


def simplex_method_tests():
    # 普通标准问题
    pb = LpProblem([-5, 5, 13, 0, 0], [[-1, 1, 3, 1, 0], [12, 4, 10, 0, 1]], [20, 90])
    s = pb.solve(SimplexMethod, show_tab=True)
    print(s)

    # 唯一解
    pb = LpProblem([2, 3, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    s = pb.solve(SimplexMethod, show_tab=True)
    print(s)

    # 无穷多解
    pb = LpProblem([2, 4, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    s = pb.solve(SimplexMethod, show_tab=True)
    print(s)

    # 无界解
    pb = LpProblem([1, 1, 0, 0], [[-2, 1, 1, 0], [1, -1, 0, 1]], [4, 2])
    s = pb.solve(SimplexMethod, show_tab=True)
    print(s)

    # 大 M 法
    pb = LpProblem([-2, -3, -1, 0, 0], [[1, 4, 2, -1, 0], [3, 2, 0, 0, -1]], [8, 6])
    s = pb.solve(SimplexMethod, show_tab=True, big_m=True)
    print(s)

    # 两阶段法
    pb = LpProblem([-2, -3, -1, 0, 0], [[1, 4, 2, -1, 0], [3, 2, 0, 0, -1]], [8, 6])
    s = pb.solve(SimplexMethod, show_tab=True, two_step=True)
    print(s)


def simplex_vectorized_tests():
    # 普通标准问题
    pb = LpProblem([-5, 5, 13, 0, 0], [[-1, 1, 3, 1, 0], [12, 4, 10, 0, 1]], [20, 90])
    s = pb.solve(SimplexVectorized)
    print(s)

    # 唯一解
    pb = LpProblem([2, 3, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    s = pb.solve(SimplexVectorized)
    print(s)

    # 无穷多解
    pb = LpProblem([2, 4, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    s = pb.solve(SimplexVectorized)
    print(s)

    # 无界解
    pb = LpProblem([1, 1, 0, 0], [[-2, 1, 1, 0], [1, -1, 0, 1]], [4, 2])
    s = pb.solve(SimplexVectorized)
    print(s)

    # 大 M 法
    pb = LpProblem([-2, -3, -1, 0, 0], [[1, 4, 2, -1, 0], [3, 2, 0, 0, -1]], [8, 6])
    s = pb.solve(SimplexVectorized)
    print(s)

    # 两阶段法
    pb = LpProblem([-2, -3, -1, 0, 0], [[1, 4, 2, -1, 0], [3, 2, 0, 0, -1]], [8, 6])
    s = pb.solve(SimplexVectorized)
    print(s)


def simplex_vectorized_debug():
    pb = LpProblem([2, 3, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    s = pb.solve(SimplexVectorized)
    print(s)


if __name__ == '__main__':
    # simplex_method_tests()
    # simplex_vectorized_debug()
    simplex_vectorized_tests()
