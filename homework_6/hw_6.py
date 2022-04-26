import random
from math import cos, sin


def first_example(x):
    return pow(x, 2) - 12 * x + 30
    # return pow(x, 2) - 2 * x + 1


def second_example(x):
    return sin(x) - cos(x)


def generate_nodes(n, x0, xn, func):
    x = [x0]
    y = [func(x0)]
    for i in range(1, n):
        aux = random.uniform(x[i - 1], xn)
        x.append(aux)
        y.append(func(aux))
    return x, y


def generate_lagrange(point, n, x, y):
    y = aitken(n, x, y)
    Ln = y[0]
    x_x = 1
    for i in range(1, n + 1):
        x_x *= (point - x[i - 1])
        Ln += y[i] * x_x
    return Ln


def aitken(n, x, y):
    for i in range(1, n + 1):
        y_aux = y.copy()
        for j in range(i, n + 1):
            y_aux[j] = (y[j] - y[j - 1]) / (x[j] - x[j - i])
        y = y_aux.copy()
    return y


print('\nf(x) = x^2 - 12 * x + 30\n')

while 1:
    try:
        n = int(input('n = '))
        if n < 3:
            raise Exception('n should be bigger that 2')
        x0 = float(input('x0 = '))
        xn = float(input('xn = '))
        point = float(input('point = '))
        if x0 < xn:
            break
        print('\tMake sure x0 < xn')
    except Exception as e:
        print(e)

func = first_example

x, y = generate_nodes(n + 1, x0, xn, func)

print('x = ', x)
print('y = ', y)
print('x0 = ', x0)
print('xn = ', xn)
print('point = ', point)

print()
Ln_point = generate_lagrange(point, n, x, y)
f_point = func(point)
print('Ln(', point, ') = ', Ln_point)
print('f(', point, ') = ', f_point)

print('| Ln( ', point, ') - f( ', point, ' ) | = ', abs(Ln_point - f_point))

print('\n ------------- \n')
print('\nf(x) = sin(x) - cos(x)\n')

while 1:
    try:
        n = int(input('n = '))
        if n < 3:
            raise Exception('n should be bigger that 2')
        x0 = float(input('x0 = '))
        xn = float(input('xn = '))
        point = float(input('point = '))
        if x0 < xn:
            break
        print('\tMake sure x0 < xn')
    except Exception as e:
        print(e)

func = second_example

x, y = generate_nodes(n + 1, x0, xn, func)

print('x = ', x)
print('y = ', y)
print('x0 = ', x0)
print('xn = ', xn)
print('point = ', point)

print()
Ln_point = generate_lagrange(point, n, x, y)
f_point = func(point)
print('Ln(', point, ') = ', Ln_point)
print('f(', point, ') = ', f_point)

print('| Ln( ', point, ') - f( ', point, ' ) | = ', abs(Ln_point - f_point))
