# Cyclotomic specialization and alternating monodromy

**Submission-stage manuscript · August 25, 2026 · 47 pages**

This repository contains the manuscript
*Cyclotomic specialization of unicritical multiplier curves and
alternating monodromy in period 4* by Dongsheng Wei.

## Main results

- For every unicritical family `f_c(z) = z^d + c`, the paper constructs a
  reduced multiplier, proves the pullback condition (PB), determines the
  maximal geometric isogeny, and gives the eventual exact cyclotomic
  factorization over the rational, maximal abelian, and fixed number-field
  bases.
- For the quadratic family, it proves the complete period-4
  Morton--Vivaldi conjecture: every delta factor `Delta_{4k,4}` is
  irreducible over the rational numbers.
- For every nontrivial root-of-unity multiplier, the associated period-4
  sextic has Galois group containing `A_6`. It remains irreducible over the
  maximal solvable extension of its cyclotomic ground field, with splitting
  group `A_6` there.
- The paper gives the exact factorization after every finite base change
  with solvable Galois closure and determines the maximal solvable subfield
  of every corresponding parabolic parameter field.
- The generic period-4 sextic has regular Galois group `S_6`, and all but
  finitely many torsion specializations retain group `S_6` over the maximal
  abelian extension of the rational numbers.

## Manuscript and certificate

- [Final PDF](final-paper/output/pdf/morton-vivaldi-period-4-and-abelian-stability.pdf)
- [LaTeX source](final-paper/main.tex)
- [Exact verification certificate](final-paper/verification/verify_all.py)
- [Outer-S5 certificate](final-paper/verification/verify_outer_s5.py)

Run the exact certificate from the repository root:

```sh
cd final-paper
python3 verification/verify_all.py
```

It uses dependency-free Python integer arithmetic to verify the companion
and exterior-power resolvent identities, three local congruences, the
period-4 parametrization, and the generic discriminant identity. It also
enumerates the exceptional `S_5` invariant and all relevant finite-group
orbits, reconstructs the outer resolvent in a rank-720 ordered splitting
algebra, verifies its discriminant from 46 exact Sylvester determinants,
and checks the unit-circle and local Taylor certificates. It uses no
floating-point or inexact numerical computation.

To build the manuscript with Tectonic:

```sh
cd final-paper
tectonic -o build main.tex
```

## Declarations

AI assistance is disclosed on the first page of the manuscript. The tools
used were OpenAI GPT-5.6 Sol Pro, OpenAI GPT-5.6 Sol Ultra through Codex,
and Anthropic Claude Opus 5 in a limited number of desktop-chat exchanges.
The author takes full responsibility for the article. The work received no
financial support, and the author declares no conflicts of interest.
