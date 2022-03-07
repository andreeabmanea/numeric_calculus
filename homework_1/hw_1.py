from cmath import pi, sin, cos, log


def find_u(base: int):
    exponent = 1
    u = pow(base, -exponent)
    while True:
        if 1.0 + u == 1.0:
            exponent -= 1
            u = pow(base, -exponent)
            break
        else:
            exponent += 1
            u = pow(base, -exponent)
    return u


# Ex. 1:
print(f"u is: {find_u(10)}\n")
print(1.0 + 1e-16 == 1.0)
print(1.0 + 1e-15 == 1.0)


def check_associativity(u: float):
    a = 1.0
    b = u / 10
    c = u / 10
    return (a + b) + c != a + (b + c)


def compute_unassoc_multiplication(u: float):
    a = 1.8
    b = u / 100
    c = u / 100
    return (a * b) * c != a * (b * c)


# Ex. 2
print(f"Sum is unassociative: {check_associativity(find_u(10))}\n")
print(
    f"Multiplication is unassociative for the given example: {compute_unassoc_multiplication(find_u(10))}\n"
)


def compute_polynom(y: float, c0: float, c1: float, c2: float, c3: float, c4: float):
    return c0 + y * (c1 + y * (c2 + y * (c3 + y * c4)))


def compute_p4_sin(y: float):
    a0 = 1805490264.690988571178600370234394843221
    a1 = 164384678.227499837726129612587952660511 * (-1)
    a2 = 3664210.647581261810227924465160827365
    a3 = 28904.140246461781357223741935980097 * (-1)
    a4 = 76.568981088717405810132543523682
    return compute_polynom(y, a0, a1, a2, a3, a4)


def compute_q4_sin(y: float):
    b0 = 2298821602.638922662086487520330827251172
    b1 = 27037050.118894436776624866648235591988
    b2 = 155791.388546947693206469423979505671
    b3 = 540.567501261284024767779280700089
    b4 = 1.0
    return compute_polynom(y, b0, b1, b2, b3, b4)


def compute_p4_cos(y: float):
    a0 = 1090157078.174871420428849017262549038606
    a1 = 321324810.993150712401352959397648541681 * (-1)
    a2 = 12787876.849523878944051885325593878177
    a3 = 150026.206045948110568310887166405972 * (-1)
    a4 = 538.333564203182661664319151379451
    return compute_polynom(y, a0, a1, a2, a3, a4)


def compute_q4_cos(y: float):
    b0 = 1090157078.174871420428867295670039506886
    b1 = 14907035.776643879767410969509628406502
    b2 = 101855.811943661368302608146695082218
    b3 = 429.772865107391823245671264489311
    b4 = 1.0
    return compute_polynom(y, b0, b1, b2, b3, b4)


def compute_p4_ln(y: float):
    a0 = 75.151856149910794642732375452928
    a1 = 134.730399688659339844586721162914 * (-1)
    a2 = 74.201101420634257326499008275515
    a3 = 12.777143401490740103758406454323 * (-1)
    a4 = 0.332579601824389206151063529971
    return compute_polynom(y, a0, a1, a2, a3, a4)


def compute_q4_ln(y: float):
    b0 = 37.575928074955397321366156007781
    b1 = 79.890509202648135695909995521310 * (-1)
    b2 = 56.215534829542094277143417404711
    b3 = 14.516971195056682948719125661717 * (-1)
    b4 = 1.0
    return compute_polynom(y, b0, b1, b2, b3, b4)


def approximate_sin(x: float):
    sin_math = sin(1 / 4 * pi * x)
    sin_approx = x * (
        compute_p4_sin(pow(x, 2))
        / (
            pow(10, -12)
            if abs(compute_q4_sin(pow(x, 2))) < pow(10, -12)
            else compute_q4_sin(pow(x, 2))
        )
    )
    print(f"argument: {x}")
    print(f"cmath sin value: {sin_math}")
    print(f"approximate sin value: {sin_approx}")
    print(f"difference in abs: {abs((sin_math - sin_approx))}")


def approximate_cos(x: float):
    cos_math = cos(1 / 4 * pi * x)
    cos_approx = compute_p4_cos(pow(x, 2)) / (
        pow(10, -12)
        if abs(compute_q4_cos(pow(x, 2))) < pow(10, -12)
        else compute_q4_cos(pow(x, 2))
    )

    print(f"argument: {x}")
    print(f"cmath cos value: {cos_math}")
    print(f"approximate cos value: {cos_approx}")
    print(f"difference in abs: {abs((cos_math - cos_approx))}")


def approximate_ln(x: float):
    z = (x - 1) / (x + 1)
    ln_math = log(x)
    ln_approx = (
        z
        * compute_p4_ln(pow(z, 2))
        / (
            pow(10, -12)
            if abs(compute_q4_ln(pow(z, 2))) < pow(10, -12)
            else compute_q4_ln(pow(z, 2))
        )
    )
    print(f"argument: {x}")
    print(f"cmath ln value: {ln_math}")
    print(f"approximate ln value: {ln_approx}")
    print(f"difference in abs: {abs((ln_math - ln_approx))}")


# Ex. 3
approximate_sin(0.5)
print()
approximate_cos(2)
print()
approximate_ln(-2)
