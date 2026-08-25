# Final manuscript

This directory contains the 47-page submission-stage manuscript
*Cyclotomic specialization of unicritical multiplier curves and the
Morton-Vivaldi conjecture in period 4* by Dongsheng Wei.

- [Final PDF](output/pdf/morton-vivaldi-period-4-and-abelian-stability.pdf)
- [LaTeX source](main.tex)
- [Exact verification certificate](verification/verify_all.py)
- [Outer-S5 certificate](verification/verify_outer_s5.py)

## Verify

```sh
python3 verification/verify_all.py
```

The dependency-free certificate uses exact Python integer arithmetic. It
checks the companion and exterior-power resolvent identities, all three
local congruences, the period-4 parametrization, and the generic
discriminant identity. It also verifies the exceptional `S_5` invariant,
the finite-group orbit table, the complete outer resolvent in a rank-720
ordered splitting algebra, its 46-point discriminant certificate, and all
local Taylor data. No floating-point or inexact numerical computation is
used.

## Build

```sh
tectonic -o build main.tex
```

## Declarations

The declarations section at the end of the manuscript gives the full AI-use
declaration, naming OpenAI GPT-5.6 Sol Pro, OpenAI GPT-5.6 Sol Ultra through
Codex, and the limited use of Anthropic Claude Opus 5 desktop chats. The
author received no financial support and declares no conflicts of interest.
