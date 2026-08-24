# Exact resolvent certificates

This directory contains a dependency-free, exact-arithmetic certificate for
the explicit `3+3` resolvent used in the abelian-stability argument.  It checks
two logically separate claims:

1. the displayed polynomial `R_X(S)` satisfies

   ```text
   det(Y I - exterior^3(C_X)) = Y^10 R_X(Y - 16/Y)
   ```

   as an identity in `Z[X,Y]`; and
2. congruences (I), (II), and (III) hold coefficientwise in `Z[u]`.

No numerical specialization is used.  The script imports only `itertools` and
`math` from the Python standard library, and all scalar arithmetic is performed
with arbitrary-precision Python integers.

## Run

From the repository root:

```sh
python3 abelian-stability/verification/verify_resolvent.py
```

The process exits with status `0` only after every check passes.  A failed
equality, nonexact division, or nonzero modular residue raises
`VerificationError` and therefore gives a nonzero exit status.

## Exact coefficient conventions

A univariate polynomial is represented by the tuple

```text
(a_0, a_1, ..., a_d)
```

for `a_0 + a_1 X + ... + a_d X^d`; trailing zero coefficients are removed.
The same representation is reused with the variable named `u` during the
congruence checks.  Thus equality means literal equality of every integer
coefficient, and reduction modulo `n` means testing every coefficient for
divisibility by `n`.

The companion matrix acts on column vectors and is

```text
        [0 0 0 0 0  16 ]
        [1 0 0 0 0   8 ]
C_X  = [0 1 0 0 0 5-X ] .
        [0 0 1 0 0  -6 ]
        [0 0 0 1 0  -4 ]
        [0 0 0 0 1  -2 ]
```

The script first independently verifies

```text
det(Y I - C_X)
 = Y^6 + 2Y^5 + 4Y^4 + 6Y^3 + (X-5)Y^2 - 8Y - 16.
```

The basis of the exterior cube is the lexicographically ordered list of the
twenty triples `I=(i<j<k)` from `{0,...,5}`.  If rows and columns are indexed by
triples `I,J`, respectively, then the constructed matrix has entry

```text
(exterior^3 C_X)_{I,J} = det((C_X)_{I,J}).
```

Every one of the 400 minors is constructed using the six-term determinant
formula over `Z[X]`.

## Determinant certificate

For any `n` by `n` matrix `A` over `Z[X]`, the script computes
`det(YI-A)` by the Faddeev--LeVerrier recurrence

```text
B_0 = I,
M_k = A B_{k-1},
c_k = -tr(M_k)/k,
B_k = M_k + c_k I       (1 <= k <= n).
```

It returns

```text
Y^n + c_1 Y^(n-1) + ... + c_n.
```

Every division by `k` is tested coefficientwise for exact divisibility in
`Z[X]`, and the terminal identity `B_n=0` is checked.  This is applied first to
`C_X` and then to the independently constructed `20` by `20` exterior-cube
matrix.

If

```text
R_X(S) = sum_{j=0}^{10} r_j(X) S^j,
```

the right-hand side of the claimed identity is independently expanded using

```text
Y^10 r_j(X)(Y-16/Y)^j
 = sum_{l=0}^j binom(j,l)(-16)^l r_j(X)
     Y^(10+j-2l).
```

The script compares the resulting coefficients of all twenty-one powers
`Y^0,...,Y^20` with the characteristic polynomial of the exterior cube.  Each
coefficient comparison is an equality in `Z[X]`, not a comparison after
specializing `X`.

The substitution `Y-16/Y` also agrees with the root interpretation: if
`rho` is a triple product, its complementary triple product is `-16/rho`, so
the paired exterior-cube eigenvalues contribute

```text
(Y-rho)(Y+16/rho) = Y^2-(rho-16/rho)Y-16.
```

## Congruence certificates

The displayed explicit coefficients `r_j(X)` are used to compute, by exact
Horner substitution in `Z[u]`, the three differences

```text
R_{u^4}((u+1)^3) - 2u^3(u+1)^27,

(partial R_X/partial X)|_{X=u^4,S=(u+1)^3} - (u+1)^26,

R_{-u^4}((u+1)^3) - 2u^3(u+1)^26.
```

The script tests every integer coefficient of these differences for
divisibility by `4`, `2`, and `4`, respectively.  It then performs the exact
division and reports the quotient degree and largest absolute quotient
coefficient as a compact reproducibility check.

## Files

- `verify_resolvent.py`: executable exact certificate.
