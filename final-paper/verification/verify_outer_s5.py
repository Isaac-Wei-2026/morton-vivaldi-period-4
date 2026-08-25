#!/usr/bin/env python3
"""Dependency-free exact certificate for the outer-S5 sextic resolvent.

The calculation takes place in the universal ordered splitting algebra of

    P_X(T) = T^6 + 2T^5 + 4T^4 + 6T^3 + (X-5)T^2 - 8T - 16.

Coefficients in Z[X] are tuples in ascending order.  An element of the
splitting algebra is a sparse dictionary on the standard 6*5*4*3*2 basis.
No numerical roots, floating-point arithmetic, factorization package, or CAS
is used.
"""

from itertools import permutations

from verify_all import (
    ONE,
    VAR,
    ZERO,
    add,
    compose as polynomial_compose,
    constant,
    evaluate_integer_polynomial,
    integer_determinant_bareiss,
    multiply,
    negate,
    poly,
    power,
    require,
    scale,
    subtract,
)


N_VARS = 5
ZERO_EXP = (0,) * N_VARS


def elem_add(left, right):
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = add(result.get(exponent, ZERO), coefficient)
        if not result[exponent]:
            del result[exponent]
    return result


def elem_negate(value):
    return {exponent: negate(coefficient) for exponent, coefficient in value.items()}


def elem_subtract(left, right):
    return elem_add(left, elem_negate(right))


def elem_scalar(value):
    return {} if not value else {ZERO_EXP: value}


class OrderedSplittingAlgebra:
    """The rank-720 ordered splitting algebra of P_X."""

    def __init__(self):
        self.relations = []
        self._monomial_cache = {}
        self.roots = []
        self._construct_roots_and_relations()

    def _clear_cache(self):
        self._monomial_cache = {}

    def reduce_monomial(self, exponent):
        """Reduce one monomial to the standard basis, with Z[X] coefficients."""
        exponent = tuple(exponent)
        cached = self._monomial_cache.get(exponent)
        if cached is not None:
            return cached

        chosen = None
        for index in range(len(self.relations) - 1, -1, -1):
            degree, _ = self.relations[index]
            if exponent[index] >= degree:
                chosen = index
                break
        if chosen is None:
            result = {exponent: ONE}
            self._monomial_cache[exponent] = result
            return result

        degree, coefficients = self.relations[chosen]
        shift = list(exponent)
        shift[chosen] -= degree
        result = {}
        # r_j^degree = -sum_{k<degree} c_k r_j^k.
        for k, coefficient_element in enumerate(coefficients):
            if not coefficient_element:
                continue
            for coefficient_exp, coefficient_poly in coefficient_element.items():
                new_exp = [
                    shift[index] + coefficient_exp[index]
                    for index in range(N_VARS)
                ]
                new_exp[chosen] += k
                reduced = self.reduce_monomial(tuple(new_exp))
                for reduced_exp, reduced_poly in reduced.items():
                    contribution = negate(multiply(coefficient_poly, reduced_poly))
                    result[reduced_exp] = add(
                        result.get(reduced_exp, ZERO), contribution
                    )
                    if not result[reduced_exp]:
                        del result[reduced_exp]
        self._monomial_cache[exponent] = result
        return result

    def multiply(self, left, right):
        result = {}
        for left_exp, left_poly in left.items():
            for right_exp, right_poly in right.items():
                product_poly = multiply(left_poly, right_poly)
                raw_exp = tuple(
                    left_exp[index] + right_exp[index]
                    for index in range(N_VARS)
                )
                reduced = self.reduce_monomial(raw_exp)
                for reduced_exp, reduced_poly in reduced.items():
                    contribution = multiply(product_poly, reduced_poly)
                    result[reduced_exp] = add(
                        result.get(reduced_exp, ZERO), contribution
                    )
                    if not result[reduced_exp]:
                        del result[reduced_exp]
        return result

    def power(self, value, exponent):
        require(exponent >= 0, "negative splitting-algebra exponent")
        result = elem_scalar(ONE)
        factor = value
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = self.multiply(result, factor)
            remaining //= 2
            if remaining:
                factor = self.multiply(factor, factor)
        return result

    def variable(self, index):
        exponent = [0] * N_VARS
        exponent[index] = 1
        return {tuple(exponent): ONE}

    def _construct_roots_and_relations(self):
        # Coefficients in ascending powers of T.
        current = [
            elem_scalar(constant(-16)),
            elem_scalar(constant(-8)),
            elem_scalar(subtract(VAR, constant(5))),
            elem_scalar(constant(6)),
            elem_scalar(constant(4)),
            elem_scalar(constant(2)),
            elem_scalar(ONE),
        ]

        for root_index in range(N_VARS):
            degree = len(current) - 1
            require(current[-1] == elem_scalar(ONE), "nonmonic quotient polynomial")
            self.relations.append((degree, current[:-1]))
            self._clear_cache()
            root = self.variable(root_index)
            self.roots.append(root)

            # Exact synthetic division by T-root.
            quotient = [{} for _ in range(degree)]
            quotient[-1] = current[-1]
            for k in range(degree - 1, 0, -1):
                quotient[k - 1] = elem_add(
                    current[k], self.multiply(root, quotient[k])
                )
            remainder = elem_add(current[0], self.multiply(root, quotient[0]))
            require(not remainder, f"nonzero synthetic remainder at root {root_index + 1}")
            current = quotient

        require(len(current) == 2 and current[1] == elem_scalar(ONE), "bad final quotient")
        sixth_root = elem_negate(current[0])
        self.roots.append(sixth_root)

        root_sum = {}
        for root in self.roots:
            root_sum = elem_add(root_sum, root)
        require(root_sum == elem_scalar(constant(-2)), "the six roots do not sum to -2")


def compose_permutations(left, right):
    """Return left after right, for image tuples."""
    return tuple(left[right[index]] for index in range(6))


def permutation_group(generators):
    identity = tuple(range(6))
    group = {identity}
    frontier = [identity]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            product = compose_permutations(generator, element)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def transposition(left, right):
    value = list(range(6))
    value[left], value[right] = value[right], value[left]
    return tuple(value)


def cycle(*entries):
    value = list(range(6))
    for left, right in zip(entries, entries[1:] + entries[:1]):
        value[left] = right
    return tuple(value)


def permute_exponent(exponent, permutation):
    result = [0] * 6
    for index, value in enumerate(exponent):
        result[permutation[index]] = value
    return tuple(result)


def outer_s5_invariant():
    # The exceptional S5=T14 and invariant from AFJR, Table 4.  The
    # permutations below are (15364), (16)(24), and (3465) in one-based
    # notation.
    generator_one = cycle(0, 4, 2, 5, 3)
    generator_two = compose_permutations(transposition(0, 5), transposition(1, 3))
    generator_three = cycle(2, 3, 5, 4)
    group = permutation_group([generator_one, generator_two, generator_three])
    require(len(group) == 120, "the proposed outer S5 does not have order 120")

    # Orbit sum of x1^2*x2^2*x3*x6 (one-based subscripts).
    base_monomial = (2, 2, 1, 0, 0, 1)
    monomials = {
        permute_exponent(base_monomial, element)
        for element in group
    }
    require(len(monomials) == 30, "outer-S5 invariant should have 30 terms")

    # Verify the full stabilizer inside S6 is exactly H.
    stabilizer = set()
    invariant = frozenset(monomials)
    for element in permutations(range(6)):
        image = frozenset(permute_exponent(term, element) for term in invariant)
        if image == invariant:
            stabilizer.add(tuple(element))
    require(stabilizer == group, "the invariant stabilizer is not exactly PGL_2(F_5)")
    return invariant, group


def evaluate_invariant(algebra, monomials, permutation):
    result = {}
    power_cache = {}
    for monomial in monomials:
        permuted = permute_exponent(monomial, permutation)
        term = elem_scalar(ONE)
        for root_index, exponent in enumerate(permuted):
            if exponent == 0:
                continue
            key = (root_index, exponent)
            if key not in power_cache:
                power_cache[key] = algebra.power(algebra.roots[root_index], exponent)
            term = algebra.multiply(term, power_cache[key])
        result = elem_add(result, term)
    return result


def explicit_a_polynomials():
    x = VAR
    return [
        subtract(x, constant(25)),
        add(add(add(power(x, 3), scale(power(x, 2), -20)), scale(x, 705)), constant(178)),
        add(add(add(add(scale(power(x, 4), 7), scale(power(x, 3), -681)), scale(power(x, 2), 16261)), scale(x, -227443)), constant(479760)),
        poly((33167184, -8149990, 1261977, -68456, 2702, -42, 1)),
        poly((2811612544, -273043148, -8878943, -1442225, 66050, 3210, -51, 3)),
        poly((-339160660992, -47208268864, 7712523489, 359206706, -56301705, 266588, -20673, 1042, 9)),
    ]


def claimed_resolvent_ascending():
    a1, a2, a3, a4, a5, a6 = explicit_a_polynomials()
    return [
        scale(a6, 16),
        scale(a5, -32),
        scale(a4, 16),
        scale(a3, 8),
        scale(a2, -8),
        scale(a1, -8),
        ONE,
    ]


def inversion_parity(permutation):
    return sum(
        permutation[i] > permutation[j]
        for i in range(6)
        for j in range(i + 1, 6)
    ) % 2


def orbit_sizes_on_cosets(subgroup, cosets, element_to_coset):
    unseen = set(range(6))
    sizes = []
    representatives = [next(iter(coset)) for coset in cosets]
    while unseen:
        point = next(iter(unseen))
        orbit = {
            element_to_coset[compose_permutations(element, representatives[point])]
            for element in subgroup
        }
        unseen -= orbit
        sizes.append(len(orbit))
    return sorted(sizes)


def verify_group_orbits(group, coset_representatives):
    cosets = [
        frozenset(compose_permutations(representative, element) for element in group)
        for representative in coset_representatives
    ]
    require(len(set(cosets)) == 6, "wrong outer-S5 coset count")
    element_to_coset = {}
    for index, coset in enumerate(cosets):
        for element in coset:
            require(element not in element_to_coset, "outer-S5 cosets overlap")
            element_to_coset[element] = index
    require(len(element_to_coset) == 720, "outer-S5 cosets do not cover S6")

    symmetric = set(permutations(range(6)))
    alternating = {
        element for element in symmetric if inversion_parity(element) == 0
    }
    exceptional_a5 = group & alternating
    pair_partition = {
        frozenset((0, 1)), frozenset((2, 3)), frozenset((4, 5))
    }
    wreath = {
        element
        for element in symmetric
        if {
            frozenset((element[left], element[right]))
            for left, right in ((0, 1), (2, 3), (4, 5))
        } == pair_partition
    }
    require(len(wreath) == 48, "wrong order for S2 wreath S3")
    expected = [
        (wreath, [2, 4], "S2 wreath S3"),
        (exceptional_a5, [1, 5], "exceptional A5"),
        (group, [1, 5], "exceptional S5"),
        (alternating, [6], "A6"),
        (symmetric, [6], "S6"),
    ]
    for subgroup, target, label in expected:
        actual = orbit_sizes_on_cosets(subgroup, cosets, element_to_coset)
        require(actual == target, f"wrong outer-coset orbit sizes for {label}")
    print(
        "PASS outer-coset group table: S2 wr S3=[2,4], "
        "exceptional A5/S5=[1,5], A6/S6=[6]"
    )


Q_POLYNOMIAL = poly((
    60081152, -12598144, 2799652, -366579, 23722, -619, 16,
))

J_POLYNOMIAL = poly((
    3812315003286784,
    -1681442076883104,
    117127525549173,
    25424297868354,
    -3439771613747,
    55056832824,
    2540918305,
    205941578,
    -3924653,
    35932,
    -1878,
    32,
))


def resolvent_coefficients_descending_at(x_value):
    a1, a2, a3, a4, a5, a6 = explicit_a_polynomials()
    return [
        1,
        -8 * evaluate_integer_polynomial(a1, x_value),
        -8 * evaluate_integer_polynomial(a2, x_value),
        8 * evaluate_integer_polynomial(a3, x_value),
        16 * evaluate_integer_polynomial(a4, x_value),
        -32 * evaluate_integer_polynomial(a5, x_value),
        16 * evaluate_integer_polynomial(a6, x_value),
    ]


def outer_resolvent_discriminant_at(x_value):
    coefficients = resolvent_coefficients_descending_at(x_value)
    derivative_coefficients = [
        (6 - index) * coefficients[index]
        for index in range(6)
    ]
    matrix = [[0] * 11 for _ in range(11)]
    for row in range(5):
        matrix[row][row:row + 7] = coefficients
    for row in range(6):
        matrix[5 + row][row:row + 6] = derivative_coefficients
    # (-1)^(6*5/2)=-1.
    return -integer_determinant_bareiss(matrix)


def verify_outer_discriminant():
    # The discriminant is weighted homogeneous of weight 30 when the
    # coefficient of Z^(6-i) has weight i.  The six X-degree/weight ratios
    # are at most 3/2, so its X-degree is at most 45.  Forty-six exact
    # specializations therefore certify the identity.
    claimed = scale(
        multiply(power(Q_POLYNOMIAL, 3), power(J_POLYNOMIAL, 2)),
        2 ** 44,
    )
    for x_value in range(46):
        require(
            outer_resolvent_discriminant_at(x_value)
            == evaluate_integer_polynomial(claimed, x_value),
            f"outer-resolvent discriminant identity fails at X={x_value}",
        )
    margin = abs(J_POLYNOMIAL[0]) - sum(
        abs(coefficient) for coefficient in J_POLYNOMIAL[1:]
    )
    require(margin == 1984823523717204, "incorrect unit-circle dominance margin")
    print(
        "PASS outer discriminant: 46 exact Sylvester determinants certify "
        "disc_Z(C_X)=2^44 Q(X)^3 J(X)^2 under the weight-45 bound"
    )
    print(f"PASS J unit-circle bound: constant-term dominance margin {margin}")


def polynomial_mod_two(value):
    return poly(coefficient % 2 for coefficient in value)


def two_adic_valuation(integer):
    require(integer != 0, "v2(0) requested in Taylor certificate")
    integer = abs(integer)
    value = 0
    while integer % 2 == 0:
        integer //= 2
        value += 1
    return value


def verify_local_data():
    x = VAR
    x_plus_one = add(x, ONE)
    expected_mod_two = [
        x_plus_one,
        multiply(x, power(x_plus_one, 2)),
        multiply(x, power(x_plus_one, 3)),
        multiply(power(x, 2), power(x_plus_one, 4)),
        multiply(power(x, 2), power(x_plus_one, 5)),
        multiply(power(x, 2), power(x_plus_one, 6)),
    ]
    expected_taylor = [
        poly((-24, 1)),
        poly((864, 668, -17, 1)),
        poly((267904, -196936, 14260, -653, 7)),
        poly((26213376, -5820800, 1072416, -58048, 2507, -36, 1)),
        poly((2528317440, -294847744, -12777920, -1146840, 81440, 2967, -30, 3)),
        poly((-378353254400, -30929592320, 8454711296, 136289280,
              -55241760, 164936, -13127, 1114, 9)),
    ]
    minus_one_units = (-13, -137, 90519, 2665647, 96165055, -4447744771)
    valuation_rows = []
    for index, a_polynomial in enumerate(explicit_a_polynomials(), start=1):
        require(
            polynomial_mod_two(a_polynomial)
            == polynomial_mod_two(expected_mod_two[index - 1]),
            f"wrong mod-2 formula for A_{index}",
        )
        translated = polynomial_compose(a_polynomial, x_plus_one)
        require(translated == expected_taylor[index - 1], f"wrong Taylor expansion for A_{index}")
        valuations = tuple(two_adic_valuation(coefficient) for coefficient in translated)
        valuation_rows.append(valuations)
        require(valuations[index] == 0, f"diagonal Taylor coefficient of A_{index} is even")
        require(
            all(valuations[j] + j > index for j in range(index)),
            f"lower Taylor term can interfere for A_{index}",
        )
        at_minus_one = evaluate_integer_polynomial(a_polynomial, -1)
        require(
            at_minus_one == (2 ** index) * minus_one_units[index - 1],
            f"wrong exact value A_{index}(-1)",
        )
    print(
        "PASS outer local data: six mod-2 identities, six exact Taylor "
        "expansions, and uniform pure-2-power valuation inequalities"
    )
    print(f"  Taylor coefficient v2 rows: {valuation_rows}")


def verify_outer_s5_resolvent():
    invariant, group = outer_s5_invariant()
    print(
        "PASS outer-S5 invariant: 30 monomials, stabilizer order 120 "
        "verified by all 720 permutations"
    )

    # AFJR's representatives: id,(56),(45),(35),(34),(46), in one-based
    # notation.
    coset_representatives = [
        tuple(range(6)),
        transposition(4, 5),
        transposition(3, 4),
        transposition(2, 4),
        transposition(2, 3),
        transposition(3, 5),
    ]
    cosets = {
        frozenset(compose_permutations(representative, element) for element in group)
        for representative in coset_representatives
    }
    require(len(cosets) == 6, "the six permutations are not distinct coset representatives")
    verify_group_orbits(group, coset_representatives)

    algebra = OrderedSplittingAlgebra()
    conjugates = [
        evaluate_invariant(algebra, invariant, representative)
        for representative in coset_representatives
    ]

    # Product in a new formal variable Z, coefficients ascending in Z.
    resolvent = [elem_scalar(ONE)]
    for conjugate in conjugates:
        updated = [{} for _ in range(len(resolvent) + 1)]
        for exponent, coefficient in enumerate(resolvent):
            updated[exponent] = elem_add(
                updated[exponent], elem_negate(algebra.multiply(coefficient, conjugate))
            )
            updated[exponent + 1] = elem_add(updated[exponent + 1], coefficient)
        resolvent = updated

    claimed = claimed_resolvent_ascending()
    for exponent, expected in enumerate(claimed):
        require(
            resolvent[exponent] == elem_scalar(expected),
            f"outer-S5 resolvent coefficient of Z^{exponent} is incorrect",
        )
    print(
        "PASS outer-S5 resolvent: all 7 coefficients reconstructed exactly "
        "in the rank-720 ordered splitting algebra"
    )

    verify_outer_discriminant()
    verify_local_data()


if __name__ == "__main__":
    verify_outer_s5_resolvent()
