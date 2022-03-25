from deepdiff import DeepDiff


def read_matrix(filename: str):

    file = open(filename, "r", closefd=True)
    dimension = int(file.readline())
    matrix = {}
    file.readline()
    new_line = file.readline().split(", ")
    while len(new_line) == 3:
        (value, line, column) = (float(new_line[0]), int(new_line[1]), int(new_line[2]))
        if not line in matrix.keys():
            matrix[line] = {}
        matrix[line].update({column: value})
        new_line = file.readline().split(", ")
    return (dimension, matrix)


def sum_matrix(file_A, file_B):
    (n_A, A) = read_matrix(file_A)
    (n_B, B) = read_matrix(file_B)
    result = {}
    if n_A != n_B:
        return False
    for i in range(0, n_A):
        for j in range(0, n_A):
            if i in A.keys() and i in B.keys():
                if j in A[i].keys() and j in B[i].keys():
                    if not i in result.keys():
                        result[i] = {}
                    sum = A[i][j] + B[i][j]
                    result[i].update({j: sum})
                elif j not in A[i].keys() and j in B[i].keys():
                    if not i in result.keys():
                        result[i] = {}
                    sum = B[i][j]
                    result[i].update({j: sum})
                elif j in A[i].keys() and j not in B[i].keys():
                    if not i in result.keys():
                        result[i] = {}
                    sum = A[i][j]
                    result[i].update({j: sum})
            elif i in A.keys() and i not in B.keys():
                result[i] = A[i]
            else:
                result[i] = B[i]

    return result


def check_sum():
    result = sum_matrix("homework_3/a.txt", "homework_3/b.txt")
    a_plus_b = read_matrix("homework_3/a_plus_b.txt")[1]
    return result == a_plus_b


def times_matrix(file_A, file_B):
    (n_A, A) = read_matrix(file_A)
    result = {}
    for i in range(0, n_A):
        if i in A.keys():
            result[i] = {}
            sum = 0
            for j in range(0, n_A):
                if j in A[i].keys():
                    print(A[i][j])
                    for k in range(0, n_A):
                        if k in A.keys():
                            sum += A[i][k] * A[k][j]
                        result[i].update({j: sum})
    return result


print(times_matrix("homework_3/a.txt", "homework_3/b.txt"))
