import json

from deepdiff import DeepDiff

epsilon = pow(10, -6)


def compare_matrix_json(C, D, n_C, n_D):

    for i in range(0, n_C):
        for j in range(0, n_C):
            if (
                str(i) in C.keys()
                and str(j) in C[str(i)].keys()
                and not (i in D.keys() and j in D[i].keys())
            ):
                return False
            if (
                i in D.keys()
                and j in D[i].keys()
                and not (str(i) in C.keys() and str(j) in C[str(i)].keys())
            ):
                return False

            if int(i) in C.keys() and int(j) in C[i].keys():
                if abs(C[i][j] - D[i][j]) > epsilon:
                    with open("difference.txt", "a") as f:
                        f.write(f"{i}, {j}, {C[i][j]}, {D[i][j]} \n")
    return True


def compare_matrix(C, D, n_C, n_D):

    for i in range(0, n_C):
        for j in range(0, n_C):
            if (
                i in C.keys()
                and j in C[i].keys()
                and not (i in D.keys() and j in D[i].keys())
            ):
                return False
            if (
                i in D.keys()
                and j in D[i].keys()
                and not (i in C.keys() and j in C[i].keys())
            ):
                return False

            if i in C.keys() and j in C[i].keys():
                if abs(C[i][j] - D[i][j]) > epsilon:
                    with open("difference.txt", "a") as f:
                        f.write(f"{i}, {j}, {C[i][j]}, {D[i][j]} \n")
    return True


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
        if column in matrix[line].keys():
            matrix[line].update({column: value + matrix[line][column]})
        else:
            matrix[line].update({column: value})
        new_line = file.readline().split(", ")
    return (dimension, matrix)


def read_json_matrix():
    with open("a_times_a.json") as json_file:
        dict = json.load(json_file)

    return dict


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
    result = sum_matrix("a.txt", "b.txt")
    a_plus_b = read_matrix("a_plus_b.txt")[1]
    return result == a_plus_b


def times_matrix(file_A, file_B):
    (n_A, A) = read_matrix(file_A)
    result = {}
    for i in range(0, n_A):
        # print(i)
        result[i] = {}
        for j in range(0, n_A):
            sum = 0
            for k in range(0, n_A):
                a = 0
                b = 0
                if i in A.keys() and k in A[i].keys():
                    a = A[i][k]
                if k in A.keys() and j in A[k].keys():
                    b = A[k][j]
                if (
                    not (i in A.keys() and k in A[i].keys())
                    and k in A.keys()
                    and i in A[k].keys()
                ):
                    a = A[k][i]
                if (
                    i in A.keys()
                    and k in A[i].keys()
                    and not (k in A.keys() and i in A[k].keys())
                ):
                    a = A[i][k]
                if (
                    not (j in A.keys() and k in A[j].keys())
                    and k in A.keys()
                    and j in A[k].keys()
                ):
                    b = A[k][j]
                if (
                    j in A.keys()
                    and k in A[j].keys()
                    and not (k in A.keys() and j in A[k].keys())
                ):
                    b = A[j][k]

                sum += a * b
            if sum != 0:
                result[i].update({j: sum})

    with open("a_times_a.json", "w") as f:
        json.dump(result, f)
    return result


def check_mul():
    result = times_matrix("a.txt", "a.txt")
    a_ori_a = read_matrix("a_ori_a.txt")[1]
    return result == a_ori_a


# print("is a+b == a_plus_b??")
# print(check_sum())
# # takes too long
# # print(check_mul())
# A = sum_matrix("a.txt", "b.txt")
# n_ab, a_plus_b = read_matrix("a_plus_b.txt")
# print("is a+b == a_plus_b??")
# print(compare_matrix(A, a_plus_b, len(A), n_ab))


print(read_matrix("test1.txt"))
print(times_matrix("test1.txt", "test1.txt"))

# print("is a_times_a = a_ori_a?")
# a_times_a = times_matrix("a.txt", "a.txt")
# a_times_a = read_json_matrix()
# # print(a_times_a)
# (n_AA, AA) = read_matrix("a_ori_a.txt")
# print(compare_matrix_json(a_times_a, AA, len(a_times_a), n_AA))
