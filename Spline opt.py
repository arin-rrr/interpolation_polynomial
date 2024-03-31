import numpy as np
import math

from matplotlib import pyplot as plt

# оба для оптимальных узлов
def f_x(x):
    return (1 / math.tan(x)) + x ** 2

def gauss1(A, b):
    for i in range(len(A)):
        if A[i][i] != 0:
            ved = A[i][i]
            for j in range(i, len(A)):
                if A[i][j] != 0:
                    A[i][j] /= ved
            if b[i] != 0:
                b[i] /= ved
            st = [0 for _ in range(len(A))]
            free_b = b[i]
            for j in range(len(A)):
                st[j] = A[i][j]
            for j in range(i+1, len(A)):
                if A[j][i] != 0:
                    d_st = st.copy()
                    free_b = b[i]
                    for k in range(len(A)):
                        d_st[k] *= A[j][i]
                    free_b *= A[j][i]
                    for k in range(len(A)):
                        A[j][k] -= d_st[k]
                    b[j] -= free_b
        else:
            for j in range(i+1, len(A)):
                if A[j][i] != 0:
                    Ai = []
                    Aj = []
                    bi = b[i]
                    bj = b[j]
                    for k in range(len(A)):
                        Ai.append(A[i][k])
                        Aj.append(A[j][k])
                    for k in range(len(A)):
                        A[i][k] = Aj[k]
                        A[j][k] = Ai[k]
                    b[i] = bj
                    b[j] = bi
                    break
            ved = A[i][i]
            for j in range(i, len(A)):
                if A[i][j] != 0:
                    A[i][j] /= ved
            if b[i] != 0:
                b[i]/= ved

            st = [0 for _ in range(len(A))]
            free_b = b[i]
            for j in range(len(A)):
                st[j] = A[i][j]
            for j in range(i+1, len(A)):
                if A[j][i] != 0:
                    d_st = st.copy()
                    free_b = b[i]
                    for k in range(len(A)):
                        d_st[k] *= A[j][i]
                    free_b *= A[j][i]
                    for k in range(len(A)):
                        A[j][k] -= d_st[k]
                    b[j] -= free_b
    x = [0 for _ in range(len(A))]
    x[-1] = b[-1]
    for i in range(len(x) - 2, -1, -1):
        res = 0
        for j in range(i+1, len(x)):
            res += A[i][j] * x[j]
        x[i] = b[i] - res
    return x

def NewtonOpt(n):
    nodes = [0 for _ in range(n+1)]
    funcInNodes = [0 for _ in range(n+1)]
    for i in range(n+1):
        nodes[i] = 0.5*(2.8 * math.cos(((2*i+1) / (2*(n+1))) * math.pi) + 3.2)
        funcInNodes[i] = f_x(nodes[i])
    razn = []
    for i in range(n+1):
        q = [0 for _ in range(n+1)]
        razn.append(q)
    for i in range(n+1):
        for j in range(n+1-i):
            if i == 0:
                razn[j][i] = funcInNodes[j]
            else:
                razn[j][i] = (razn[j + 1][i - 1] - razn[j][i - 1]) / (nodes[i + j] - nodes[j])

    res = []
    for i in range(n+1):
        if i == 0:
            sv = razn[i][i]
            q = [[]]
            res.append([sv, q])
        else:
            sv = razn[0][i]
            q = []
            for j in range(i):
                q.append([1, -nodes[j]])
            res.append([sv, q])
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

def S32opt(n):
    nodes = [0 for _ in range(n+1)]
    funcInNodes = [0 for _ in range(n+1)]
    for i in range(n+1):
        nodes[i] = 0.5*(2.8 * math.cos(((2*i+1) / (2*(n+1))) * math.pi) + 3.2)
    nodes.sort()
    for i in range(n+1):
        funcInNodes[i] = f_x(nodes[i])

    h = [0 for _ in range(n)]
    for i in range(n):
        h[i] = nodes[i+1] - nodes[i]
    H = []
    for i in range(n-1):
        a = [0 for _ in range(n-1)]
        H.append(a)
    for i in range(n-1):
        if i == 0:
            H[i][i] = 2* (h[i] + h[i+1])
            H[i][i+1] = h[i+1]
        else:
            if i == n-2:
                H[i][i] = 2*(h[i] + h[i+1])
                H[i][i-1] = h[i]
            else:
                H[i][i - 1] = h[i]
                H[i][i] = 2 * (h[i] + h[i + 1])
                H[i][i + 1] = h[i + 1]

    gamma = [0 for _ in range(n-1)]
    for i in range(1, n):
        gamma[i-1] = 6 * ((funcInNodes[i + 1] - funcInNodes[i]) / h[i] - (funcInNodes[i] - funcInNodes[i - 1]) / h[i - 1])
    vt_pr = gauss1(H, gamma)
    second_der = [0 for _ in range(n+1)]
    for i in range(n+1):
        if i == 0 or i == n:
            second_der[i] = 0;
        else:
            second_der[i] = vt_pr[i-1]
    first_der = [0 for _ in range(n)]
    for i in range(n):
        first_der[i] = (funcInNodes[i + 1] - funcInNodes[i]) / h[i] - (second_der[i + 1] * h[i]) / 6 - (second_der[i] * h[i]) / 3

    koef = []
    for i in range(n):
        a = [0 for _ in range(4)]
        a[0] = (second_der[i + 1] - second_der[i]) / (6 * h[i])
        a[1] = second_der[i] / 2
        a[2] = first_der[i]
        a[3] = funcInNodes[i]
        koef.append(a)
    return koef

def S32(nums, x):
    n = len(nums)
    nodes = [0 for _ in range(n+1)]
    for i in range(n+1):
        nodes[i] = 0.5*(2.8 * math.cos(((2*i+1) / (2*(n+1))) * math.pi) + 3.2)
    nodes.sort()
    if x <= nodes[0]:
        y = 0
    elif x >= nodes[-1]:
        y = -1
    else:
        for i in range(n):
            if x >= nodes[i] and x <= nodes[i + 1]:
                y = i
                break
    res = nums[y][0]* ((x-nodes[y]) ** 3) + nums[y][1]*((x- nodes[y]) ** 2) + nums[y][2] * (x-nodes[y]) + nums[y][3]
    return res

k = int(input('Количество узлов '))

x = np.arange(0.2, 3, 0.01)
y_newton = []
y_spl = []

res1 = NewtonOpt(k)
res2 = S32opt(k)

for i in x:
    y_newton.append(abs(newton(res1, i) - f_x(i)))
    y_spl.append(abs(S32(res2, i) - f_x(i)))
print(y_newton)
print(y_spl)


plt.plot(x, y_newton, label=rf'N(x)')
plt.plot(x, y_spl, label=rf'S32(x)')
plt.xlabel('Ось х')  # Подпись для оси х
plt.ylabel('Ось y')  # Подпись для оси y
plt.title(f'Распределение абсолютной погрешности ньютона и S32 (оптимальные узлы) с {k} узлами')  # Название
plt.legend(loc='best', fontsize=10)
plt.grid()
plt.show()