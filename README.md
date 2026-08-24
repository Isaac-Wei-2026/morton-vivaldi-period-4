# Period-4 Morton--Vivaldi Irreducibility and Abelian Stability

**Version 2.1 · manuscript snapshot · August 25, 2026**

This repository contains a proof of the complete period-4 case of the
Morton--Vivaldi irreducibility conjecture for the quadratic family

\[
f_c(z)=z^2+c.
\]

## Main results

The main algebraic theorem states that, for every integer `k >= 2`,

\[
\Gamma_k(T)=(-1)^{\varphi(k)}T^{2\varphi(k)}\Phi_k(g(T))
\]

is irreducible over the rational numbers, where

\[
g(T)=-T^4-2T^3-4T^2-6T+5+8/T+16/T^2.
\]

Using the normalization of the period-4 multiplier curve, the paper deduces
the irreducibility of every delta factor `Delta_{4k,4}`.

Version 2.1 also proves the abelian-stability strengthening: for
every nontrivial root of unity `zeta`, the sextic fiber `F_zeta` remains
irreducible over the maximal abelian extension of `Q(zeta)`.

## Version 2.1 manuscript

- [Read the versioned PDF](final-paper/output/pdf/morton-vivaldi-period-4-and-abelian-stability.pdf)
- [LaTeX source](final-paper/main.tex)
- [Exact verification certificate](final-paper/verification/verify_all.py)

Earlier standalone manuscripts are retained in [paper/](paper/) and
[abelian-stability/](abelian-stability/).

## Verify and build

With Python 3 and Tectonic installed, run from the repository root:

```sh
cd final-paper
python3 verification/verify_all.py
tectonic -X compile main.tex --outdir build
```

The compiled file is written to `final-paper/build/main.pdf`.

The core irreducibility proof was obtained using OpenAI Sol Pro.
