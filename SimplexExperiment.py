from linprog.problem import LpProblem
from linprog.simplex_method import SimplexMethod
from linprog.simplex_vectorized import SimplexVectorized
from linprog.dual_simplex import DualSimplex


def simplex_method_tests():
    """
    单纯形表法测试
    """
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
    """
    矩阵表示的单纯形法测试
    """
    # 上课的例题 for debug
    pb = LpProblem([2, 3, 0, 0, 0], [[1, 2, 1, 0, 0], [4, 0, 0, 1, 0], [0, 4, 0, 0, 1]], [8, 16, 12])
    s = pb.solve(SimplexVectorized)
    print(s)

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


def dual_simplex_tests():
    """
    对偶单纯形法测试
    """
    pb = LpProblem([-9, -12, -15, 0, 0, 0], [[-2, -2, -1, 1, 0, 0], [-2, -3, -1, 0, 1, 0], [-1, -1, -5, 0, 0, 1]], [-10, -12, -14])
    s = pb.solve(DualSimplex)
    print(s)

    pb = LpProblem([-5, 5, 13, 0, 0], [[-1, 1, 3, 1, 0], [16, 0, -2, -4, 1]], [30, -30])
    s = pb.solve(DualSimplex)
    print(s)

    pb = LpProblem([-5, 5, 13, 0, 0, 0], [[-1, 1, 3, 1, 0, 0], [16, 0, -2, -4, 1, 0], [5, 0, -4, -3, 0, 1]], [20, 10, -10])
    s = pb.solve(DualSimplex)
    print(s)


def experiment():
    pb = LpProblem([3, -1, -1, 0, 0], [[1, -2, 1, 1, 0], [-4, 1, 2, 0, -1], [-2, 0, 1, 0, 0]], [11, 3, 1])

    print("\n\n>> 单纯形表法，打印单纯形表，使用两阶段法: ")
    s1 = pb.solve(SimplexMethod, show_tab=True, two_step=True)
    print(s1)

    print("\n\n>> 矩阵表示单纯形，默认使用大M法: ")
    s2 = pb.solve(SimplexVectorized)
    print(s2)

    print("\n\n>> 对偶单纯形法: ")
    pbd = LpProblem([-9, -12, -15, 0, 0, 0], [[-2, -2, -1, 1, 0, 0], [-2, -3, -1, 0, 1, 0], [-1, -1, -5, 0, 0, 1]], [-10, -12, -14])
    s3 = pbd.solve(DualSimplex)
    print(s3)


if __name__ == '__main__':
    # simplex_method_tests()
    # simplex_vectorized_tests()
    # dual_simplex_tests()
    experiment()
