#!/usr/bin/env python3
"""Exact certificates for the period-4 3+3 resolvent.

This program uses only Python integers.  A polynomial is a tuple of its
integer coefficients in ascending order.  In particular, no floating-point
arithmetic, numerical specialization, or computer-algebra package is used.
"""

from itertools import combinations
from math import comb


class VerificationError(RuntimeError):
    """Raised when an exact certificate check fails."""


def require(condition, message):
    if not condition:
        raise VerificationError(message)


# Univariate polynomials over Z, stored in ascending coefficient order.
ZERO = ()
ONE = (1,)
VAR = (0, 1)


def poly(values):
    values = list(values)
    while values and values[-1] == 0:
        values.pop()
    return tuple(values)


def constant(value):
    return ZERO if value == 0 else (int(value),)


def add(left, right):
    size = max(len(left), len(right))
    out = [0] * size
    for index in range(size):
        if index < len(left):
            out[index] += left[index]
        if index < len(right):
            out[index] += right[index]
    return poly(out)


def negate(value):
    return tuple(-coefficient for coefficient in value)


def subtract(left, right):
    return add(left, negate(right))


def scale(value, scalar):
    if scalar == 0 or not value:
        return ZERO
    return tuple(scalar * coefficient for coefficient in value)


def multiply(left, right):
    if not left or not right:
        return ZERO
    out = [0] * (len(left) + len(right) - 1)
    for i, a_i in enumerate(left):
        if a_i == 0:
            continue
        for j, b_j in enumerate(right):
            if b_j:
                out[i + j] += a_i * b_j
    return poly(out)


def power(value, exponent):
    require(exponent >= 0, "negative polynomial exponent")
    result = ONE
    factor = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, factor)
        remaining //= 2
        if remaining:
            factor = multiply(factor, factor)
    return result


def derivative(value):
    return poly(index * value[index] for index in range(1, len(value)))


def compose(outer, inner):
    """Return outer(inner), exactly, by Horner evaluation."""
    result = ZERO
    for coefficient in reversed(outer):
        result = add(multiply(result, inner), constant(coefficient))
    return result


def divide_by_integer_exact(value, divisor):
    require(divisor != 0, "division by zero")
    quotients = []
    for coefficient in value:
        quotient, remainder = divmod(coefficient, divisor)
        require(
            remainder == 0,
            f"coefficient {coefficient} is not divisible by {divisor}",
        )
        quotients.append(quotient)
    return poly(quotients)


def degree(value):
    return len(value) - 1


def identity_matrix(size):
    matrix = [[ZERO for _ in range(size)] for _ in range(size)]
    for index in range(size):
        matrix[index][index] = ONE
    return matrix


def matrix_multiply(left, right):
    """Exact matrix product over Z[X], with zero-entry skipping."""
    row_count = len(left)
    inner_count = len(right)
    column_count = len(right[0])
    require(all(len(row) == inner_count for row in left), "matrix size mismatch")
    result = [[ZERO for _ in range(column_count)] for _ in range(row_count)]
    for i in range(row_count):
        for k, left_entry in enumerate(left[i]):
            if not left_entry:
                continue
            for j, right_entry in enumerate(right[k]):
                if right_entry:
                    result[i][j] = add(
                        result[i][j], multiply(left_entry, right_entry)
                    )
    return result


def matrix_trace(matrix):
    result = ZERO
    for index in range(len(matrix)):
        result = add(result, matrix[index][index])
    return result


def characteristic_polynomial(matrix):
    """Return det(YI-A), descending in Y, by Faddeev--LeVerrier.

    If the output is [1,c_1,...,c_n], it represents
    Y^n+c_1 Y^(n-1)+...+c_n.  Every division by k is checked exactly in Z[X].
    """
    size = len(matrix)
    require(size > 0, "empty matrix")
    require(all(len(row) == size for row in matrix), "matrix is not square")
    b_matrix = identity_matrix(size)
    coefficients = [ONE]
    for k in range(1, size + 1):
        product = matrix_multiply(matrix, b_matrix)
        c_k = divide_by_integer_exact(negate(matrix_trace(product)), k)
        coefficients.append(c_k)
        for index in range(size):
            product[index][index] = add(product[index][index], c_k)
        b_matrix = product
    require(
        all(not entry for row in b_matrix for entry in row),
        "Faddeev--LeVerrier terminal matrix is nonzero",
    )
    return coefficients


def companion_matrix():
    """Companion matrix of F_X in the column-vector convention."""
    matrix = [[ZERO for _ in range(6)] for _ in range(6)]
    for column in range(5):
        matrix[column + 1][column] = ONE
    # Negatives of the coefficients from T^0 through T^5.
    last_column = [
        constant(16),
        constant(8),
        subtract(constant(5), VAR),
        constant(-6),
        constant(-4),
        constant(-2),
    ]
    for row, value in enumerate(last_column):
        matrix[row][5] = value
    return matrix


def determinant_3_by_3(matrix):
    positive = add(
        add(
            multiply(multiply(matrix[0][0], matrix[1][1]), matrix[2][2]),
            multiply(multiply(matrix[0][1], matrix[1][2]), matrix[2][0]),
        ),
        multiply(multiply(matrix[0][2], matrix[1][0]), matrix[2][1]),
    )
    negative = add(
        add(
            multiply(multiply(matrix[0][2], matrix[1][1]), matrix[2][0]),
            multiply(multiply(matrix[0][1], matrix[1][0]), matrix[2][2]),
        ),
        multiply(multiply(matrix[0][0], matrix[1][2]), matrix[2][1]),
    )
    return subtract(positive, negative)


def exterior_cube(matrix):
    """Matrix of exterior^3(matrix) on lexicographically ordered triples."""
    require(len(matrix) == 6 and all(len(row) == 6 for row in matrix), "not 6x6")
    basis = list(combinations(range(6), 3))
    result = [[ZERO for _ in basis] for _ in basis]
    for out_index, rows in enumerate(basis):
        for in_index, columns in enumerate(basis):
            minor = [[matrix[row][column] for column in columns] for row in rows]
            result[out_index][in_index] = determinant_3_by_3(minor)
    return basis, result


def explicit_resolvent():
    """Coefficients r_j(X) of R_X(S)=sum_{j=0}^{10} r_j(X) S^j."""
    x = VAR
    x2 = power(x, 2)
    x3 = power(x, 3)
    result = [ZERO] * 11
    result[10] = ONE
    result[9] = constant(6)
    result[8] = scale(add(x, constant(35)), 4)
    result[7] = scale(add(add(x2, scale(x, -6)), constant(389)), 2)
    result[6] = add(
        add(add(x3, scale(x2, -15)), scale(x, 667)), constant(5683)
    )
    result[5] = scale(add(add(x2, scale(x, -2)), constant(337)), 96)
    result[4] = scale(
        add(add(add(x3, scale(x2, -14)), scale(x, 477)), constant(800)),
        64,
    )
    result[3] = scale(add(add(x2, scale(x, -21)), constant(-368)), -1536)
    result[2] = scale(add(add(x2, scale(x, 35)), constant(-48)), 8192)
    result[1] = scale(multiply(subtract(x, constant(16)), add(x, constant(4))), -65536)
    result[0] = scale(power(subtract(x, constant(16)), 2), -65536)
    return result


def y10_resolvent_substitution(resolvent):
    """Coefficients in Y of Y^10 R_X(Y-16/Y), ascending in Y."""
    result = [ZERO] * 21
    for s_exponent, coefficient in enumerate(resolvent):
        # (Y-16/Y)^j = sum_l binom(j,l)(-16)^l Y^(j-2l).
        for l in range(s_exponent + 1):
            y_exponent = 10 + s_exponent - 2 * l
            multiplier = comb(s_exponent, l) * ((-16) ** l)
            result[y_exponent] = add(
                result[y_exponent], scale(coefficient, multiplier)
            )
    return result


def resolvent_evaluate(resolvent, x_value, s_value):
    """Evaluate R_X(S) at polynomial values X=x_value, S=s_value."""
    result = ZERO
    for coefficient in reversed(resolvent):
        result = add(multiply(result, s_value), compose(coefficient, x_value))
    return result


def verify_coefficientwise_divisibility(value, modulus, label):
    require(modulus > 0, "nonpositive modulus")
    for exponent, coefficient in enumerate(value):
        require(
            coefficient % modulus == 0,
            f"{label}: coefficient of u^{exponent} is {coefficient}, not 0 mod {modulus}",
        )
    quotient = divide_by_integer_exact(value, modulus)
    inspected = max(1, len(value))
    max_quotient = max((abs(coefficient) for coefficient in quotient), default=0)
    print(
        f"PASS {label}: {inspected} coefficients checked mod {modulus}; "
        f"quotient degree {degree(quotient)}, max |coefficient| {max_quotient}"
    )


def verify_companion_and_resolvent():
    companion = companion_matrix()
    expected_f = [
        ONE,
        constant(2),
        constant(4),
        constant(6),
        subtract(VAR, constant(5)),
        constant(-8),
        constant(-16),
    ]
    actual_f = characteristic_polynomial(companion)
    require(actual_f == expected_f, "the 6x6 matrix is not a companion matrix of F_X")
    print("PASS companion: det(YI-C_X)=F_X(Y), all 7 Y-coefficients checked")

    basis, wedge = exterior_cube(companion)
    require(len(basis) == 20 and len(wedge) == 20, "exterior cube has wrong size")
    nonzero_entries = sum(bool(entry) for row in wedge for entry in row)
    print(
        "PASS exterior cube: constructed 20x20 matrix from all 3x3 minors "
        f"({nonzero_entries} nonzero entries)"
    )

    q_descending = characteristic_polynomial(wedge)
    q_ascending = list(reversed(q_descending))
    claimed = y10_resolvent_substitution(explicit_resolvent())
    require(len(q_ascending) == 21 and len(claimed) == 21, "wrong Y-degree")
    for exponent in range(21):
        require(
            q_ascending[exponent] == claimed[exponent],
            f"resolvent identity fails at the coefficient of Y^{exponent}",
        )
    print(
        "PASS resolvent identity: det(YI-exterior^3 C_X)="
        "Y^10 R_X(Y-16/Y), all 21 Y-coefficients equal in Z[X]"
    )


def verify_three_congruences():
    resolvent = explicit_resolvent()
    u = VAR
    u4 = power(u, 4)
    u_plus_one = add(u, ONE)
    s_value = power(u_plus_one, 3)

    # (I)
    left_i = resolvent_evaluate(resolvent, u4, s_value)
    right_i = scale(multiply(power(u, 3), power(u_plus_one, 27)), 2)
    verify_coefficientwise_divisibility(
        subtract(left_i, right_i), 4, "congruence (I)"
    )

    # (II): differentiate every coefficient r_j(X), then evaluate.
    derivative_resolvent = [derivative(coefficient) for coefficient in resolvent]
    left_ii = resolvent_evaluate(derivative_resolvent, u4, s_value)
    right_ii = power(u_plus_one, 26)
    verify_coefficientwise_divisibility(
        subtract(left_ii, right_ii), 2, "congruence (II)"
    )

    # (III)
    left_iii = resolvent_evaluate(resolvent, negate(u4), s_value)
    right_iii = scale(multiply(power(u, 3), power(u_plus_one, 26)), 2)
    verify_coefficientwise_divisibility(
        subtract(left_iii, right_iii), 4, "congruence (III)"
    )


def main():
    print("Exact arithmetic domain: Z[X] and Z[u] (Python arbitrary-precision integers)")
    verify_companion_and_resolvent()
    verify_three_congruences()
    print("ALL EXACT CERTIFICATE CHECKS PASSED")


if __name__ == "__main__":
    main()
