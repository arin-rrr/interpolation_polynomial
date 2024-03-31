import numpy as np
import math

from matplotlib import pyplot as plt

k = int(input('Количество узлов '))


def f_x(x):
    return x ** 2 + 1 / math.tan(x)


def NewtonRavn(n):
    nodes = [0 for _ in range(n + 1)]
    funcInNodes = [0 for _ in range(n + 1)]
    for i in range(n + 1):
        nodes[i] = 0.2 + i * (2.8 / n)
        funcInNodes[i] = f_x(nodes[i])

    razn = []
    for i in range(n + 1):
        razn.append([0 for _ in range(n + 1)])

    for i in range(n + 1):
        for j in range(n + 1 - i):
            if i == 0:
                razn[j][i] = funcInNodes[j]
            else:
                razn[j][i] = (razn[j + 1][i - 1] - razn[j][i - 1]) / (nodes[i + j] - nodes[j])

    res = []
    for i in range(n + 1):
        if i == 0:
            res.append([razn[i][i], [[]]])
        else:
            q = []
            for j in range(i):
                q.append([1, -nodes[j]])
            res.append([razn[0][i], q])
    return res


def NewtonOpt(n):
    nodes = [0 for _ in range(n + 1)]
    funcInNodes = [0 for _ in range(n + 1)]
    for i in range(n + 1):
        nodes[i] = 0.5 * (2.8 * math.cos(math.pi * (2 * i + 1) / (2 * (n + 1))) + 3.2)
        funcInNodes[i] = f_x(nodes[i])
    razn = []
    for i in range(n + 1):
        razn.append([0 for _ in range(n + 1)])

    for i in range(n + 1):
        for j in range(n + 1 - i):
            if i == 0:
                razn[j][i] = funcInNodes[j]
            else:
                razn[j][i] = (razn[j + 1][i - 1] - razn[j][i - 1]) / (nodes[i + j] - nodes[j])

    res = []
    for i in range(n + 1):
        if i == 0:
            res.append([razn[i][i], [[]]])
        else:
            q = []
            for j in range(i):
                q.append([1, -nodes[j]])
            res.append([razn[0][i], q])
    return res

def newton(nums, x):
    result = 0
    for j in range(len(nums)):
        if j == 0:
            result+=nums[j][0]
        else:
            pr_res = nums[j][0]
            for k in range(len(nums[j][1])):
                pr_res *= (x + nums[j][1][k][1])
            result += pr_res
    return result

res1 = NewtonRavn(k)
res2 = NewtonOpt(k)

x = np.arange(0.2, 3, 0.01)
y_ravn = []
y_opt = []

for i in x:
    y_ravn.append(newton(res1, i))
    y_opt.append(newton(res2, i))

plt.plot(x, x ** 2 + 1 / np.tan(x), label=r'f(x)')
plt.plot(x, y_ravn, label=rf'N{k}(x)')
plt.plot(x, y_opt, label=rf'N_opt{k}(x)')
plt.xlabel('Ось х')  # Подпись для оси х
plt.ylabel('Ось y')  # Подпись для оси y
plt.title(f'Ньютон с {k} узлами интерполяции')  # Название
plt.legend(loc='best', fontsize=10)
plt.grid()
plt.show()