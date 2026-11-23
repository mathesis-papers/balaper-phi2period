# Polynomial-flow extensions of Su's theorem on Bala's period conjectures

Haobo Ma — ChronoAI Pte Ltd, Singapore  
Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore

License: Creative Commons Attribution 4.0 International (CC BY 4.0).

## Mathematical contents

The note gives an alternative finite-jet proof of Bala's five eventual totient-square conjectures, already proved by Lien-Hung Su, *Quadratic Derivations and Uniform Modular Periods for Six Inverse-Function Sequences*, version 1, 31 October 2026, Theorem 3.2 and Corollary 3.3, DOI `10.5281/zenodo.23835577`.

For the integral quadratic derivations `D_a(s)=c^a`, `D_a(c)=s*c^(a-1)`, `c^2-s^2=1`, the extension covers every integer exponent `a>=1` and every integral polynomial observable at even derivative indices. Modulo `m>=3`, the shift `phi(m)^2/2` is eventual. The explicit onset `J(m)` is the maximum of `ceil(p^(2r-1)/2)` over odd prime powers exactly dividing `m` and `2^(r-1)` for its 2-primary factor. Modulo 2 the column period is one from index 1.

The written corollaries cover a finite-difference description of residue rows, A281182 interleaving, eight neighboring inverse-function sequences, the affirmative answer to the secant-operator question in A281436, the factorial normalization in A281428, the nested-derivative offsets in A281427 and A281435, the common inverse representations and odd-exponent algebraic identities, and the derivative-polynomial triangle A394012 with its A391180 second-column formula and an eventual period. Su's stronger A281428 vanishing and totient-period theorem is cited separately with its precise scope.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a standard TeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is included in the TeX source; BibTeX is unnecessary.

## Exact calculations

Python 3.10 or later and its standard library suffice:

```sh
python3 verify.py --output reproduced-checks.json
cmp checks.json reproduced-checks.json
```

The deterministic output includes:

- 21 coefficients in each of the five primary families, first-eight-term OEIS comparisons, and a separate ODE coefficient calculation through derivative index 40.
- Seven least-period certificates, with finite matrices, zero shift differences and smaller-period witnesses. The Cayley–Hamilton reduction in the note proves their all-index implications.
- Five geometric-row certificates, the all-index identity A281181(n)=1 modulo 12, and the constant tail A281184(n)=3 modulo 9 for n>=1.
- 60 polynomial-operator identities, and 108 quadratic-operator cases, each checking both the period identity and the finite-difference annihilator for exponents 1 through 9 at twelve prime powers.
- Thirteen coefficients in each of eight neighboring families, with their first eight OEIS terms; exact secant-operator comparisons through derivative index 12; four factorial-normalization comparisons through derivative index 12; and twelve complete rows of the triangle, compared with its source recurrence and second-column identity.
- Two exact algebraic elimination identities, the quotient-series identity for T_4 through degree 13, and finite second-column period comparisons at fourteen moduli from 2 through 15.

The finite operator cases complement the all-parameter written proofs. The pinned Lean declarations linked individually in the PDF cover the original five totient-square assertions and the arbitrary-polynomial prime-power lemma. The all-exponent extension, improved onset, residue-row extension, sibling identities and numerical certificate reduction have written proofs in the note. `metadata.json` records the corresponding scopes.

## Sources

- Su: https://doi.org/10.5281/zenodo.23835577
- Five posing entries: https://oeis.org/A281180, https://oeis.org/A281181, https://oeis.org/A281183, https://oeis.org/A281184, https://oeis.org/A281440
- Bala's inverse-function identity: https://oeis.org/A281181/a281181.pdf
- Dominici's nested derivatives: https://doi.org/10.1155/S0161171203303291
- Danesh's weighted Stirling-transform theorem: https://doi.org/10.33774/coe-2026-4d6jr
- Additional formulas and triangle: https://oeis.org/A281436, https://oeis.org/A281428, https://oeis.org/A281427, https://oeis.org/A281435, https://oeis.org/A394012, https://oeis.org/A391180

The source archive contains the exposition, calculation code, deterministic results, bibliographic metadata and a SHA-256 manifest. Third-party papers and database snapshots are available from the cited sources.
