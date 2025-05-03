import numpy as np
import matplotlib.pyplot as plt

# 全局变量
M = 50  # x方向等分
N = 50  # y方向等分


# 定义SOR方法计算的函数
def SOR(psi, omega, dx, dy, w, tolerance):
    beta = dx / dy
    maxe = 1
    iter = 0
    while (maxe > tolerance):
        iter += 1
        maxe = 0
        for i in range(1, M - 1):
            for j in range(1, N - 1):
                psi1 = psi[i][j]
                psi[i][j] = w / 2 / (1 + beta ** 2) * (
                            psi[i + 1][j] + psi[i - 1][j] + beta ** 2 * (psi[i][j + 1] + psi[i][j - 1]) + dx ** 2 *
                            omega[i][j]) + (1 - w) * psi[i][j]
                maxe = max(maxe, abs(psi[i][j] - psi1))
    return iter  # 返回迭代次数


# 定义找最佳松弛因子的函数
def SOR_w_opt(p, omega, dx, dy, w, tolerance):
    min_i = SOR(p.copy(), omega, dx, dy, w, tolerance)
    w_opt = w
    for w1 in np.arange(w + 0.01, 1.86, 0.01):
        i = SOR(p.copy(), omega, dx, dy, w1, tolerance)
        if (min_i > i):
            min_i = i
            w_opt = w1
    return w_opt


# 初始化参数
dx = 1.0 / M
dy = 1.0 / N
dt = 0.001  # 时间
nu = 0.001  # 运动粘度
max_iter = 2000  # 迭代时间步数
tolerance = 1e-7  # SOR迭代的收敛限
w = 1.86  # 初始松弛因子
w_opt = w

# 初始化数组
psi = np.zeros((M, N))  # 流函数
omega = np.zeros((M, N))  # 涡量
u = np.zeros((M, N))
v = np.zeros((M, N))

# 初始时刻边界的涡量，Thom公式
for i in range(1, M - 1):
    omega[i][N - 1] = -2 * np.sin(np.pi * i * dx) ** 2 / dy

omega_new = omega.copy()  # 辅助计算

# 初始时刻速度场
for i in range(0, M):
    u[i][N - 1] = np.sin(np.pi * i * dx) ** 2

# 设置实时绘图
plt.ion()
fig, ax = plt.subplots(figsize=(8, 6))
x = np.linspace(0, 1, M)
y = np.linspace(0, 1, N)
X, Y = np.meshgrid(x, y, indexing='ij')  # 生成网格坐标


for t in range(max_iter):
    # 计算n+1时刻内点的涡量
    for i in range(1, M - 1):
        for j in range(1, N - 1):
            # 计算扩散项，五点离散laplace算子
            diffusion = (omega[i + 1, j] + omega[i - 1, j] - 2 * omega[i, j]) / dx ** 2 + (
                        omega[i, j + 1] + omega[i, j - 1] - 2 * omega[i, j]) / dy ** 2

            # 计算对流项
            convection = u[i, j] * (omega[i + 1, j] - omega[i - 1, j]) / (2 * dx) + v[i, j] * (
                        omega[i, j + 1] - omega[i, j - 1]) / (2 * dy)

            # 更新内点涡量
            omega_new[i][j] = omega[i][j] + dt * (-convection + nu * diffusion)

    omega = omega_new.copy()  # 更新涡量

    # 计算w_opt
    if (t == 0):
        print("正在计算最佳松弛因子")
        w_opt = SOR_w_opt(psi, omega, dx, dy, w, tolerance)
        print("最佳松弛因子:", w_opt)

    # 更新流函数，用SOR方法解泊松方程
    _ = SOR(psi, omega, dx, dy, w_opt, tolerance)

    # 更新速度场
    for i in range(1, M - 1):
        for j in range(1, N - 1):
            u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * dy)
            v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * dx)

    # 更新上下边界上的涡量
    for i in range(0, M):
        omega[i][0] = -2 * (psi[i][1] - psi[i][0]) / dy ** 2
        omega[i][N - 1] = -2 * (psi[i][N - 2] - psi[i][N - 1] + np.sin(np.pi * i * dx) ** 2 * dy) / dy ** 2
    # 更新左右边界上的涡量
    for j in range(0, N):
        omega[0][j] = -2 * (psi[1][j] - psi[0][j]) / dx ** 2
        omega[M - 1][j] = -2 * (psi[M - 2][j] - psi[M - 1][j]) / dx ** 2

    # 绘制流线图
    if t % 10 == 0:
        ax.cla()
        levels_ = np.arange(psi.min(), psi.max()+0.001, 0.001)
        contour = ax.contour(X, Y, psi, 
                        levels=levels_, 
                        colors='k', 
                        linestyles='solid', 
                        linewidths=0.5)
        ax.set_title(f'Time Step: {t}') 
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        plt.pause(0.01)