# SimplexLinprog
单纯形算法实现。

## 项目结构

```
.
├── LICENSE
├── README.md
├── SimplexExperiment.py		算法测试 
└── linprog
    ├── __init__.py
    ├── dual_simplex.py			对偶单纯形算法实现 DualSimplex
    ├── problem.py			线性规划问题类 LpProblem
    ├── simplex_method.py		单纯形表法实现 SimplexMethod
    ├── simplex_vectorized.py		单纯形算法矩阵表示实现 SimplexVectorized
    ├── solve.py			线性规划问题的解类 LpSolve
    └── solver.py			线性规划问题解法抽象类 LpSolver
```

## 接口说明

在开始之前，您需要了解该项目提供的三个基本模型类：

- `LpProblem`：表示一个**标准型**线性规划问题：

  ```
  max z = c^T * x
  s.t. a*x = b, x >= 0, b > 0
  ```

  【注】使用 `DualSimplex` 算法求解的时候，还有注意保证问题可用对偶单纯形法。

- `LpSolve`：表示一个线性规划问题的解，包含以下内容：

  ```
  success: 最优化成功 : True or False
  description: 解的描述 : 唯一最优解/无穷多最优解/无界解/无可行解... 
  solve: 最优解 : [x1, x2, ...] or None
  target: 最优目标函数值 : z
  ```

- `LpSolver`：线性规划问题的解法，具体实现有三种：

  ```
  SimplexMethod		# 单纯形表法
  SimplexVectorized	# 基于矩阵表示的单纯形法
  DualSimplex		# 对偶单纯形法
  ```

**具体使用**：

1. 初始化一个线性规划标准型问题（`LpProblem`）

   ```python
   lpProblem = LpProblem(c, a, b)
   """
   :param c: 目标函数系数
   :param a: 系数矩阵
   :param b: 右端常数
   """
   ```

2. 指定具体算法，求解问题

   ```python
   solve = lpProblem.solve(solver, **kwargs)
   """
   :param solver: LpSolver 的具体实现，直接传类名即可，三选一：
   	SimplexMethod		: 单纯形表法
   	SimplexVectorized	: 基于矩阵表示的单纯形法
   	DualSimplex		: 对偶单纯形法
   :param kwargs: 可选参数
   	base_idx=[...] (default []):    指定初始基，缺省则由算法自行确定
   	show_tab=True  (default False): 打印单纯形表，仅对 solver=SimplexMethod 适用
   	two_step=True  (default False): 无单位阵作初始基时，使用两阶段法
   	big_m=True     (default True):  无单位阵作初始基时，使用大 M 法
   :return: 问题的解 LpSolve
   """
   ```

3. 查看解

   ```python
   print(solve)
   # LpSolve 实现了对__str__的重载，可以打印出易读的解
   # 或者，可以调取 LpSolve 对象的具体属性：
   solve.success     # 是否得到了最优解
   solve.description # 解的描述
   solve.solve       # 最优解向量
   solve.target      # 最优目标函数值
   ```

E.g.

```python
# 导入需要的命名
from linprog.problem import LpProblem

from linprog.simplex_method import SimplexMethod
# from linprog.simplex_vectorized import SimplexVectorized
# from linprog.dual_simplex import DualSimplex

# 初始化一个线性规划标准型问题
pb = LpProblem([-5, 5, 13, 0, 0], [[-1, 1, 3, 1, 0], [12, 4, 10, 0, 1]], [20, 90])
# 用单纯形表法求解，打印单纯形表
s = pb.solve(SimplexMethod, show_tab=True)
# 打印解
print(s)
```

运行结果：

```
   x_3	  20.0	  -1.0	   1.0	[3.00]	   1.0	     0	6.6666	
   x_4	  90.0	  12.0	   4.0	  10.0	     0	   1.0	   9.0	
      	      	  -5.0	   5.0	  13.0	     0	     0	      	
----------------
   x_2	6.6666	-0.333	[0.33]	   1.0	0.3333	     0	20.000	
   x_4	23.333	15.333	0.6666	     0	-3.333	   1.0	34.999	
      	      	-0.666	0.6666	     0	-4.333	     0	      	
----------------
   x_1	20.000	  -1.0	   1.0	   3.0	   1.0	     0	      	
   x_4	9.9999	  16.0	     0	-2.000	  -4.0	   1.0	      	
      	      	     0	     0	-2.000	  -5.0	     0	      	
----------------

最优化成功	: True
解的描述	: 无穷多最优解
最优解	: [ 0. 20.  0.  0. 10.]
最优目标函数值	: 100.00000000000001
```

打印出的单纯形表格式如下：

|      |      |       |       |
| :--: | :--: | :---: | :---: |
| x_B  |  b   |   A   | theta |
|      |      | sigma |       |

## TODO

- [ ] 重构，抽象 SimplexMethod、SimplexVectorized、DualSimplex，避免代码重复
- [ ] 异常输入检测、处理
- [ ] 输入任意线性规划问题，自动化为标准型
- [ ] CLI、GUI 交互界面
- [ ] 模仿 scipy.optimize.linprog，重写一套更好的实现

## 开放源代码

MIT License

Copyright (c) 2020 CDFMLR
