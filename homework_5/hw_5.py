from numpy import cumsum
import numpy as np
from math import sqrt

k_max = 20000
EPSILON = pow(10, -6)


def getI_N(n):
    I_N = []
    for i in range(0, n):
        for j in range(0, n):
            if i == j:
                I_N.append(1)
            else:
                I_N.append(0)

    return np.array(I_N).reshape(n, n)


def get_position_in_array(i, j, N):
    return int(i * N - (i - 1) * i / 2 + j - i)


def convert_matrix_to_array(matrix):
    array = []
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i >= j:
                array.append(matrix[i][j])

    return array


def get_diagonal_indices(N: int):
    list = [x for x in range(1, N + 1)]
    return cumsum(list)


def get_row_col(index: int, diagonal_indices: list):
    for diagonal_index in range(len(diagonal_indices)):
        if diagonal_indices[diagonal_index] > index:
            row = diagonal_index
            break
    if row > 0:
        col = index - diagonal_indices[row - 1]
    else:
        col = 0
    return (row, col)


def get_p_and_q(array, N):
    diagonal_indices = get_diagonal_indices(N)
    maxim_value = -100000
    maxim_index = -1
    for i in range(len(array)):
        (row, col) = get_row_col(i, diagonal_indices)
        if row == col:
            continue
        if abs(array[i]) > maxim_value:
            maxim_index = i
            maxim_value = abs(array[i])

    return get_row_col(maxim_index, diagonal_indices)


def compute_teta(array, p, q, N):
    pq = get_position_in_array(p, q, N)
    pp = get_position_in_array(p, p, N)
    qq = get_position_in_array(q, q, N)
    alfa = (array[pp] - array[qq]) / (2 * array[pq])  # cotg(2teta)
    if alfa >= 0:
        t = -alfa + sqrt(pow(alfa, 2) + 1)
    else:
        t = -alfa - sqrt(pow(alfa, 2) + 1)

    c = 1 / (sqrt(1 + pow(t, 2)))
    s = t / (sqrt(1 + pow(t, 2)))
    return (c, s, t)


def is_diagonal_matrix(array, N):
    non_diagonal_elements = []
    diagonal_indices = get_diagonal_indices(N)
    for i in range(len(array)):
        (row, col) = get_row_col(i, diagonal_indices)
        if row == col:
            continue
        non_diagonal_elements.append(array[i])
    return all(v <= EPSILON for v in non_diagonal_elements)


def compute_rotation(N, p, q, c, s):
    R = np.zeros((N, N))
    for i in range(len(N)):
        for j in range(len(N)):
            if i == j and i != p and i != q:
                R[i][j] = 1
            elif i == j and (i == p or i == q):
                R[i][j] = c
            elif i == p and j == q:
                R[i][j] = s
            elif i == q and j == p:
                R[i][j] = -s
            else:
                R[i][j] = 0
    return R


def compute_new_A(array, p, q, c, s, t, N):
    new_A = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i >= j:
                new_A[i][j] = array[get_position_in_array(i, j, N)]
                new_A[j][i] = new_A[i][j]
    for j in range(N):
        if j != p and j != q:
            new_A[p][j] = (
                c * array[get_position_in_array(p, j, N)]
                + s * array[get_position_in_array(q, j, N)]
            )
            new_A[j][p] = (
                c * array[get_position_in_array(p, j, N)]
                + s * array[get_position_in_array(q, j, N)]
            )
            new_A[q][j] = (
                -s * array[get_position_in_array(p, j, N)]
                + c * array[get_position_in_array(q, j, N)]
            )
            new_A[j][q] = (
                -s * array[get_position_in_array(p, j, N)]
                + c * array[get_position_in_array(q, j, N)]
            )

    new_A[p][p] = (
        pow(c, 2) * array[get_position_in_array(p, p, N)]
        + pow(s, 2) * array[get_position_in_array(q, q, N)]
        + 2 * c * s * array[get_position_in_array(p, q, N)]
    )
    new_A[q][q] = (
        pow(s, 2) * array[get_position_in_array(p, p, N)]
        + pow(c, 2) * array[get_position_in_array(q, q, N)]
        - 2 * c * s * array[get_position_in_array(p, q, N)]
    )
    new_A[p][q] = (pow(c, 2) - pow(s, 2)) * array[
        get_position_in_array(p, q, N)
    ] + c * s * (
        array[get_position_in_array(q, q, N)] - array[get_position_in_array(p, p, N)]
    )
    new_A[q][p] = (pow(c, 2) - pow(s, 2)) * array[
        get_position_in_array(p, q, N)
    ] + c * s * (
        array[get_position_in_array(q, q, N)] - array[get_position_in_array(p, p, N)]
    )
    return new_A


def compute_new_U(U, p, q, s, c, N):
    new_U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i >= j:
                new_U[i][j] = U[i][j]
                new_U[j][i] = new_U[i][j]
    for i in range(N):
        new_U[i][p] = c * U[i][p] + s * U[i][q]
        new_U[i][q] = -s * U[i][p] + c * U[i][q]
    return new_U


def jacobi_algorithm(array, N):
    k = 0
    U = getI_N(N)
    (p, q) = get_p_and_q(array, N)
    (c, s, t) = compute_teta(array, p, q, N)
    while not is_diagonal_matrix(array, N) and k <= k_max:

        array = convert_matrix_to_array(compute_new_A(array, p, q, c, s, t, N))
        U = compute_new_U(U, p, q, s, c, N)
        (p, q) = get_p_and_q(array, N)
        (c, s, t) = compute_teta(array, p, q, N)
        k = k + 1
    return (array, U)


# A = [[0, 1, 3, 6], [1, 2, 4, 7], [3, 4, 5, 8], [6, 7, 8, 9]]
A = [[1, 1, 2], [1, 1, 2], [2, 2, 2]]
# A = [[0, 0, 1], [0, 0, 1], [1, 1, 1]]
# A = [[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
copy_A = A.copy()
array = convert_matrix_to_array(A)
N = len(A)
arr, U = jacobi_algorithm(array, N)
final_A = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if i >= j:
            final_A[i][j] = arr[get_position_in_array(i, j, N)]
            final_A[j][i] = final_A[i][j]
prod = np.dot(np.dot(U.transpose(), copy_A), U)
print(prod)
print(final_A.tolist())
