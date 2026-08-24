# Version 2.1 manuscript snapshot

This directory contains the merged manuscript proving the complete
period-4 case of the Morton--Vivaldi irreducibility conjecture and its
abelian-stability strengthening.

- [Versioned PDF](output/pdf/morton-vivaldi-period-4-and-abelian-stability.pdf)
- [LaTeX source](main.tex)
- [Exact verification certificate](verification/verify_all.py)

## Verify

From this directory, run:

```sh
python3 verification/verify_all.py
```

The exact certificate uses only Python integer arithmetic.  It verifies
the exterior-power resolvent identity, the three local polynomial
congruences, and the period-4 multiplier-curve parametrization identity.

## Build

```sh
tectonic -X compile main.tex --outdir build
```

The compiled file is written to `build/main.pdf`.
