import numpy as np

MAX_STEP = 10000
EPSILON = pow(10, -6)


def read_matrix(filename: str):
    # am adaugat vectorul diagonal care retine elementele de pe diagonala principala
    file = open(filename, "r", closefd=True)
    dimension = int(file.readline())
    matrix = {}
    file.readline()
    new_line = file.readline().split(", ")
    diagonal = []
    while len(new_line) == 3:
        (value, line, column) = (float(new_line[0]), int(new_line[1]), int(new_line[2]))
        if line == column:
            diagonal.append(value)
        if not (line, column) in matrix.keys():
            matrix[line, column] = value
        else:
            matrix[line, column] += matrix[line, column]
        new_line = file.readline().split(" , ")
    return (dimension, matrix, diagonal)


def read_vector_b_column(filename: str):
    file = open(filename, "r", closefd=True)
    lines = file.readlines()
    b = [float(line.strip()) for line in lines]
    return b


def read_and_check_input(filename_A, filename_b):
    (n_A, A, diagonal) = read_matrix(filename_A)
    b = read_vector_b_column(filename_b)
    for element in diagonal:
        if element < pow(10, -6):
            raise Exception("There are some 0 elements in the main diagonal")
    return (n_A, A, b)


def solve_with_jacobi(A, b, n_A):
    x_c = [0.0] * n_A
    x_p = [0.0] * n_A
    convergent = True
    step = 0
    delta = pow(10, -4)
    while delta >= EPSILON and delta <= pow(10, 8) and step <= MAX_STEP:
        x_p = x_c.copy()
        for position, _ in A.items():
            step += 1
            x_c[position[0]] = 0.0 * x_c[position[0]]
            first_sum = 0
            second_sum = 0
            if position[1] < position[0]:
                first_sum += A[position] * x_c[position[1]]
            if position[1] >= (position[0] + 1):
                second_sum += A[position] * x_p[position[1]]
            numerator = b[position[0]] - (first_sum + second_sum)
            x_c[position[0]] += (1.0 * numerator) / A[position[0], position[0]]
            diff = [0.0] * n_A
            for i in range(n_A):
                diff[i] = x_c[i] - x_p[i]
            delta = np.linalg.norm(diff)
    if delta < EPSILON:
        print(f"Solution: {x_c} was found after {step} iterations")
        return x_c
    else:
        convergent = False
        return "Divergent"


def check_solution(A, n_A, b, x_sol):
    result = [0.0] * n_A
    difference = [0.0] * n_A
    for position, _ in A.items():
        result[position[0]] += A[position] * x_sol[position[1]]
    for i in range(n_A):
        difference[i] = result[i] - b[i]
    return f"The norm for the solution is: {np.linalg.norm(difference)}"


# takes too long
# (n_A, A, b) = read_and_check_input("homework_4/a_1.txt", "homework_4/b_1.txt")
# this example was taken from: https://colab.research.google.com/drive/1buZaoFdC1NRAISoZI7o0wYP4TXtRR595
# the result from there are the same with ours
(n_A, A, b) = read_and_check_input("homework_4/test_a.txt", "homework_4/test_b.txt")
solution = solve_with_jacobi(A, b, n_A)
print(check_solution(A, n_A, b, solution))
