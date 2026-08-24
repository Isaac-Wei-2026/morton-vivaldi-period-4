# The Period-4 Morton--Vivaldi Irreducibility Problem

This repository contains a proof of the complete period-4 case of the
Morton--Vivaldi irreducibility conjecture for the quadratic family

\[
f_c(z)=z^2+c.
\]

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

Version 2.0 adds a standalone strengthening: for every nontrivial root of
unity `zeta`, the sextic fiber `F_zeta` remains irreducible over the maximal
abelian extension of `Q(zeta)`.

## Files

- `paper/main.tex` - LaTeX source of the paper.
- `paper/references.bib` - bibliography database.
- `output/pdf/morton-vivaldi-period-4.pdf` - compiled paper.
- `abelian-stability/draft.tex` - source of the Version 2.0 strengthening.
- `abelian-stability/output/pdf/abelian-stability-v2.pdf` - compiled Version
  2.0 paper.
- `abelian-stability/verification/verify_resolvent.py` - exact certificate for
  the resolvent identity and local congruences.

## Local compilation

With Tectonic installed, run from the repository root:

```sh
cd paper
tectonic --outdir ../output/pdf main.tex
mv ../output/pdf/main.pdf ../output/pdf/morton-vivaldi-period-4.pdf
```

The core irreducibility proof was obtained using OpenAI Sol Pro.
