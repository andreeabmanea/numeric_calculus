import random
from sympy import Symbol, sympify

K_MAX = 10000

EPSILON = pow(10, -6)


def f_x(function_str: str, value: float):
    x = Symbol("x")
    sympified = sympify(function_str)
    return float(sympified.subs(x, value))


def steffensen_method(F: str, G_method: int, h: float, x_init: float):
    x = x_init
    k = 0
    while True:
        g_x = compute_g(G_method, F, x, h)
        x_2 = x + g_x
        g_x_2 = compute_g(G_method, F, x_2, h)
        if abs(g_x_2 - g_x) < EPSILON:
            return x
        delta_x = pow(g_x, 2) / (g_x_2 - g_x)
        x = x - delta_x
        k = k + 1
        if delta_x < EPSILON or k > K_MAX or delta_x >= pow(10, 8):
            break
    if compute_double_derivative(F, x, h) > 0:
        print("x* is minimum point")
    else:
        print("x* is not minimum point")
    if delta_x < EPSILON:
        return x
    else:
        return "divergent"


def G1(F, x, h):
    F_x = f_x(F, x)
    F_x_h = f_x(F, (x - h))
    F_x_2h = f_x(F, (x - 2 * h))
    return (3 * F_x - 4 * F_x_h + F_x_2h) / (2 * h)


def G2(F, x, h):
    F_x_plus_2h = f_x(F, (x + 2 * h))
    F_x_plus_h = f_x(F, (x + h))
    F_x_minus_h = f_x(F, (x - h))
    F_x_minus_2h = f_x(F, (x - 2 * h))
    return (-F_x_plus_2h + 8 * F_x_plus_h - 8 * F_x_minus_h + F_x_minus_2h) / (12 * h)


def compute_g(method: int, F, x, h):
    if method == 1:
        return G1(F, x, h)
    else:
        return G2(F, x, h)


def compute_double_derivative(F, x, h):
    F_x_plus_2h = f_x(F, (x + 2 * h))
    F_x_plus_h = f_x(F, (x + h))
    F_x_minus_h = f_x(F, (x - h))
    F_x_minus_2h = f_x(F, (x - 2 * h))
    F_x = f_x(F, x)
    return (
        -F_x_plus_2h + 16 * F_x_plus_h - 30 * F_x + 16 * F_x_minus_h - F_x_minus_2h
    ) / (12 * pow(h, 2))


print("F(x) = x^2 + sin(x)")
result = steffensen_method("x**2 + sin(x)", 1, pow(10, -6), -0.5)
print(f"G1: {result}")
result = steffensen_method("x**2 + sin(x)", 2, pow(10, -6), -0.5)
print(f"G2: {result}")
print()
print("F(x) = 1/3 * (x ** 3) - 2 * (x ** 2) + 2 * x + 3")
result = steffensen_method(
    "1/3 * (x ** 3) - 2 * (x ** 2) + 2 * x + 3", 1, pow(10, -6), random.uniform(3, 3.5)
)
print(f"G1: {result}")
result = steffensen_method(
    "1/3 * (x ** 3) - 2 * (x ** 2) + 2 * x + 3", 2, pow(10, -6), random.uniform(3, 3.5)
)
print(f"G2: {result}")
print()
print("F(x) =  x ** 4 - 6 * (x ** 3) + 13 * (x ** 2) - 12 * x + 4")
result = steffensen_method(
    "x ** 4 - 6 * (x ** 3) + 13 * (x ** 2) - (12 * x) + 4",
    1,
    pow(10, -6),
    random.uniform(1, 2),
)
print(f"G1: {result}")
result = steffensen_method(
    "x ** 4 - 6 * (x ** 3) + 13 * (x ** 2) - (12 * x) + 4",
    2,
    pow(10, -6),
    random.uniform(1, 2),
)
print(f"G2: {result}")
