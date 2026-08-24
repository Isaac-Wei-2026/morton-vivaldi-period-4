# Abelian stability of the period-4 fibers

This directory contains the standalone Version 2.0 manuscript proving the
abelian-stable strengthening of the period-4 Morton--Vivaldi fiber theorem.
It is separate from the companion manuscript in `paper/`.

## Main theorem

For every nontrivial root of unity `zeta`, with `K = Q(zeta)`, the sextic
`F_zeta` is irreducible over the maximal abelian extension `K^ab`.
Equivalently, for a root `alpha`,

```text
K(alpha) intersection K^ab = K.
```

The manuscript also proves the resulting factorization theorem over finite
abelian base fields and over `Q^ab`, transfers it to the period-4 delta
factors, determines the maximal abelian subfield of the corresponding
parameter fields, and identifies the unique exceptional torsion fiber.

## Files

- `draft.tex`: integrated standalone LaTeX manuscript.
- `drafts/core_lemmas.tex`: Newton polygons, the local quartic factor, and the
  abstract sextic criterion.
- `drafts/resolvent_local.tex`: explicit 3+3 resolvent and local no-root proof.
- `drafts/corollaries.tex`: abelian-base, delta-factor, parameter-field, and
  torsion-Hilbert consequences.
- `references.bib`: bibliography.
- `verification/verify_resolvent.py`: dependency-free exact certificate for
  the resolvent identity and congruences (I)--(III).
- `output/pdf/abelian-stability-v2.pdf`: final visually verified PDF.

## Reproduce

From this directory:

```sh
python3 verification/verify_resolvent.py
tectonic -X compile draft.tex --outdir build
```

The certificate uses exact integer polynomial arithmetic and exits nonzero on
any failed equality, nonexact division, or failed coefficientwise congruence.
The LaTeX manuscript compiles without warnings.
