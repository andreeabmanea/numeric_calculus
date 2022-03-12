import numpy as np


def example_A():
    return np.array([[2, 0, 1], [0, 2, 1], [4, 4, 6]])


def example_B():
    return [[5], [1], [14]]


def generate_A(n: int):
    return np.random.randint(10, size=(n, n))


def generate_B(n: int):
    return np.split(np.random.randint(50, size=n), n)


def is_singular(matrix: np.array):
    if np.linalg.det(matrix) == 0:
        return True
    return False


def choose_pivot(matrix: np.array, l: int):
    pivot_candidates = []
    for row in matrix:
        pivot_candidates.append(row[l - 1])
    pivot_candidates = [abs(x) for x in pivot_candidates]

    return int(pivot_candidates.index(np.amax(pivot_candidates)) + 1)


def swap_lines(matrix: np.array, i: int, l: int):
    matrix[[i - 1, l]] = matrix[[l, i - 1]]
    return matrix


def add_B_to_A(A: np.array, B: np.array):
    return np.append(A, B, axis=1)


def gauss_elimination(A, B, n, eps):
    l = 1
    A = add_B_to_A(A, B)
    pivot = choose_pivot(A, l)
    A = swap_lines(A, pivot, l)
    while l < n - 1 and abs(A[l - 1][l - 1]) > eps:
        for i in range(l, n):
            f = A[i][l - 1] / A[l - 1][l - 1]
            for j in range(l, n + 1):
                A[i][j] = A[i][j] - f * A[l - 1][j]
                A[i][l - 1] = 0
        l += 1
        pivot = choose_pivot(A, l)
        A = swap_lines(A, pivot, l)
    if abs(A[l - 1][l - 1]) <= eps:
        print("Singular matrix")
        return
    else:
        det = 1
        for i in range(0, n):
            det *= A[i][i]
        if det == 0:
            print("A is singular")
        else:
            B = A[:, n]
            A = A[:, 0:n]
            print(solve_upper_triangular_matrix(A, B, n))


def solve_upper_triangular_matrix(A, b, n):
    x_gauss = np.zeros(n)

    for i in range(n - 1, -1, -1):
        temp = b[i]
        for j in range(n - 1, i, -1):
            temp -= x_gauss[j] * A[i, j]

        x_gauss[i] = temp / A[i, i]
    return x_gauss


gauss_elimination(example_A(), example_B(), 3, pow(10, -5))
gauss_elimination(generate_A(5), generate_B(5), 5, pow(10, -5))
