import numpy as np
import math

from matplotlib import pyplot as plt

k = int(input('Количество узлов '))
def f_x(x):
    return x ** 2 + 1 / math.tan(x)

def Lagrange(nums, x):
    result = 0
    for j in range(len(nums)):
        pr_res = 1
        for k in range(len(nums[j][1])):
            pr_res *= x + nums[j][1][k][1]
        pr_res *= nums[j][0]
        result += pr_res
    return result

def LagrRavn(n):
    nodes = [0 for _ in range(n+1)]
    funcInNodes = [0 for _ in range(n+1)]
    for i in range(n+1):
        nodes[i] = 0.2+i*(2.8 / n)
        funcInNodes[i] = f_x(nodes[i])

    res = []
    for i in range(n+1):
        znam = 1
        work = []
        for j in range(n+1):
            if i!= j:
                znam *= (nodes[i] - nodes[j])
                work.append([1, -nodes[j]])
        res.append([funcInNodes[i]/znam, work])
    return res

def LagrOpt(n):
    nodes = [0 for _ in range(n+1)]
    funcInNodes = [0 for _ in range(n+1)]

    for i in range(n+1):
        nodes[i] = 0.5*(2.8 * math.cos(math.pi * (2*i+1) / (2*(n+1))) + 3.2)
        funcInNodes[i] = f_x(nodes[i])

    res = []
    for i in range(k):
        znam = 1
        work = []
        for j in range(n+1):
            if i != j:
                znam *= (nodes[i] - nodes[j])
                work.append([1, -nodes[j]])
        mn = f_x(nodes[i]) / znam
        res.append([mn, work])
    return res

res1 = LagrRavn(k)
res2 = LagrOpt(k)

x = np.arange(0.2, 3, 0.01)
y_ravn = []
y_opt = []

for i in x:
    y_ravn.append(Lagrange(res1, i))
    y_opt.append(Lagrange(res2, i))

plt.plot(x, x ** 2 + 1 / np.tan(x), label=r'f(x)')
plt.plot(x, y_ravn, label=rf'L{k}(x)')
plt.plot(x, y_opt, label=rf'L_opt{k}(x)')
plt.xlabel('Ось х')  # Подпись для оси х
plt.ylabel('Ось y')  # Подпись для оси y
plt.title(f'Лагранж с {k} узлами интерполяции')  # Название
plt.legend(loc='best', fontsize=10)
plt.grid()
plt.show()