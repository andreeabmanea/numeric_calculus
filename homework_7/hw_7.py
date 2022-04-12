import random

from numpy import number

EPSILON = pow(10, -6)
K_MAX = 10000


def compute_r(coefficients: list):
    coefficients = [abs(x) for x in coefficients]
    A = max(coefficients)
    return (coefficients[0] + A) / coefficients[0]


def evaluate_poly_with_horner(x, a):
    result = 0
    for i in range(len(a) - 1, -1, -1):
        result = a[i] + (x * result)
    return int(result)


def dehghan_method(coefficients):
    R = compute_r(coefficients)
    x_0 = random.uniform(-R, R)
    k = 0
    while True:
        if abs(evaluate_poly_with_horner(x_0, coefficients)) <= EPSILON / 10:
            delta_x = 0
        else:
            poly_x_k = evaluate_poly_with_horner(x_0, coefficients)
            poly_x_k_plus_poly = evaluate_poly_with_horner(x_0 + poly_x_k, coefficients)
            poly_x_k_minus_poly = evaluate_poly_with_horner(
                x_0 - poly_x_k, coefficients
            )
            y_k = x_0 - (2 * pow(poly_x_k, 2)) / (
                poly_x_k_plus_poly - poly_x_k_minus_poly
            )
            poly_y_k = evaluate_poly_with_horner(y_k, coefficients)
            delta_x = (2 * poly_x_k * (poly_x_k + poly_y_k)) / (
                poly_x_k_plus_poly - poly_x_k_minus_poly
            )
        x_0 = x_0 - delta_x
        k += 1
        if delta_x <= EPSILON or k <= K_MAX or delta_x > pow(10, 8):
            break
    if delta_x < EPSILON:
        return x_0
    else:
        return "divergenta"


def create_output(number_of_runs, coefficients):
    results = []
    for i in range(number_of_runs):
        x = dehghan_method(coefficients)
        new = True
        for j in range(len(results)):
            if abs(results[j] - x) < EPSILON:
                new = False
        if new:
            results.append(x)

    print(results)
    output_file = open("homework_7/output.txt", "w")

    for result in results:
        output_file.write(str(result) + "\n")

    output_file.close()


# coefficients: a0 (termenul liber), a1, ... an
create_output(100, [-6, 11, -6, 1])
