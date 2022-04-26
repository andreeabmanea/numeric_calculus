from numpy import cumsum
import numpy as np
from math import sqrt

np.set_printoptions(suppress=True)

k_max = 1000
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


def get_p_and_q(array):
    maxim_value = -100000
    maxim_index_i = -1
    maxim_index_j = -1
    for ii in range(len(array)):
        for jj in range(len(array[ii])):
            if ii < jj:
                if abs(array[ii][jj]) > maxim_value:
                    maxim_index_i = ii
                    maxim_index_j = jj
                    maxim_value = abs(array[ii][jj])

    return maxim_index_i, maxim_index_j


def compute_teta(array, p, q, N):
    # pq = get_position_in_array(p, q, N)
    # pp = get_position_in_array(p, p, N)
    # qq = get_position_in_array(q, q, N)
    alfa = (array[p][p] - array[q][q]) / (2 * array[p][q])  # cotg(2teta)
    if alfa >= 0:
        t = -alfa + sqrt(pow(alfa, 2) + 1)
    else:
        t = -alfa - sqrt(pow(alfa, 2) + 1)

    c = 1 / (sqrt(1 + pow(t, 2)))
    s = t / (sqrt(1 + pow(t, 2)))
    return c, s, t


def is_diagonal_matrix(array):
    for i in range(len(array)):
        for j in range(len(array[i])):
            if i == j:
                continue
            if array[i][j] > EPSILON:
                return False
    return True


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
                new_A[i][j] = array[i][j]
                new_A[j][i] = new_A[i][j]
    for j in range(N):
        if j != p and j != q:
            new_A[p][j] = (
                    c * array[p][j]
                    + s * array[q][j]
            )
            new_A[j][p] = (
                    c * array[p][j]
                    + s * array[q][j]
            )
            new_A[q][j] = (
                    -s * array[p][j]
                    + c * array[q][j]
            )
            new_A[j][q] = (
                    -s * array[p][j]
                    + c * array[q][j]
            )

    new_A[p][p] = (
            pow(c, 2) * array[p][p]
            + pow(s, 2) * array[q][q]
            + 2 * c * s * array[p][q]
    )
    new_A[q][q] = (pow(s, 2) * array[p][p] + pow(c, 2) * array[q][q] - 2 * c * s * array[p][q])
    new_A[p][q] = (pow(c, 2) - pow(s, 2)) * array[p][q] + c * s * (array[q][q] - array[p][p])
    new_A[q][p] = (pow(c, 2) - pow(s, 2)) * array[p][q] + c * s * (array[q][q] - array[p][p])
    return new_A


def compute_new_U(U, p, q, s, c, N):
    new_U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            new_U[i][j] = U[i][j]
    for i in range(N):
        new_U[i][p] = c * U[i][p] + s * U[i][q]
        new_U[i][q] = -s * U[i][p] + c * U[i][q]
    return new_U


def compute_S(singular_values, rank, P, N):
    new_S = np.zeros((N, P))
    for i in range(rank):
        if singular_values[i] > 0:
            new_S[i][i] = 1 / singular_values[i]

    return new_S


def jacobi_algorithm(array, N):
    global P
    if P != N:
        return '\n !----- p != N -----! \n'
    k = 0
    U = getI_N(N)
    (p, q) = get_p_and_q(array)
    (c, s, t) = compute_teta(array, p, q, N)
    while not is_diagonal_matrix(array) and k <= k_max:
        array = compute_new_A(array, p, q, c, s, t, N)
        U = compute_new_U(U, p, q, s, c, N)
        (p, q) = get_p_and_q(array)
        (c, s, t) = compute_teta(array, p, q, N)
        k = k + 1
    return array, U


def get_eigen_values(A):
    eigen_values = []
    for i in range(len(A)):
        for j in range(len(A[i])):
            if i == j:
                eigen_values.append(A[i][j])
    return eigen_values


def compute_least_squares_pseudo_inverse(A):
    return np.linalg.inv(A.transpose() @ A) @ A.transpose()


def compute_Moore_Penrose_pseudo_inverse(A):
    U, S, V = np.linalg.svd(A)
    S = compute_S(S, np.linalg.matrix_rank(A), len(A), len(A[0]))
    return V.T @ S @ U.T


# A = [[0, 1, 3, 6], [1, 2, 4, 7], [3, 4, 5, 8], [6, 7, 8, 9]]
# A = [[1, 1, 2], [1, 1, 2], [2, 2, 2]]
# A = [[0, 0, 1], [0, 0, 1], [1, 1, 1]]

A = [[0.68482643, 0.21107947, 0.1649274, -0.04874136, 1.187186],
     [0.21107947, 0.55784448, 0.07125934, 0.18908876, 0.27793555],
     [0.1649274, 0.07125934, 0.88050656, 0.11329253, 0.67456082],
     [-0.04874136, 0.18908876, 0.11329253, 0.55033413, -0.29052412],
     [1.187186, 0.27793555, 0.67456082, -0.29052412, 5.14630661]]
# A = [[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
A = np.array(A)
copy_A = A.copy()
P = len(A) + 1
N = len(A[0])

if P == N:
    arr, U = jacobi_algorithm(A, N)

    print("\nA final:")
    print(arr)

    print("\nEigen values:")
    eigen_values = get_eigen_values(arr)
    print(eigen_values)

    print("\nU: ")
    print(U)

    print("\n||A init * U - U ADiag ||")
    print(np.linalg.norm(copy_A @ U - U @ arr))

    lib_eigen = np.linalg.eigh(copy_A)
    print("\nEigen values from NP library:")
    print(lib_eigen[0])

    print("\nEigen vector from NP library:")
    print(lib_eigen[1])

    sum = 0
    for i in range(len(eigen_values)):
        min = 99999
        for k in range(len(lib_eigen[0])):
            if abs(eigen_values[i] - lib_eigen[0][k]) < min:
                min = abs(eigen_values[i] - lib_eigen[0][k])
        sum += min

    print("\nThe sum difference between eigen_values and eigen_values_from_NP is: ")
    print(sum)
else:
    P -= 1
    A = [[0., 0., 1.],
         [0., 0., 1.],
         [1., 1., 1.],
         [1., 1., 1.]]
    A = np.array(A)
    svd = np.linalg.svd(A)
    print('\nSingular values of matrix A using SVD:')
    print(svd[1])

    print('\nRank matrix for A using SVD:')
    rank = 0
    for i in svd[1]:
        if i > EPSILON:
            rank += 1
    print(rank)

    print('\nRank matrix for A using SVD lib:')
    rank = np.linalg.matrix_rank(A)
    print(rank)

    min_plus = 9999
    max_plus = -9999
    for i in svd[1]:
        if i > pow(10, -6):
            if min_plus > i:
                min_plus = i
        if max_plus < i:
            max_plus = i
    print('\nConditional number for matrix A:')
    print(max_plus / min_plus)

    A_i = compute_Moore_Penrose_pseudo_inverse(A)
    print('\nMoore-Penrose pseudo inverse of A is: ')
    print(A_i)

    print("\nMoore-Penrose pseudo inverse of A from NP:")
    print(np.linalg.pinv(A))

    A_j = compute_least_squares_pseudo_inverse(A)
    print('\nLeast square pseudo inverse of A is: ')
    print(A_j)

    print('\n|| A_i - A_j || 1')
    print(np.linalg.norm(np.subtract(A_i, A_j), ord=1))
