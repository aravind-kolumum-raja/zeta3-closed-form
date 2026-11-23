# 100-Search Literature Archive on $\zeta(3)$ (Apéry's Constant)

Comprehensive index of the **102 targeted web research queries** (10 thematic clusters × 10 queries + 2 surgical formula verification queries) executed on November 24, 2026, spanning 1734 through November 2026.

---

## Cluster 1 (Searches 1–10): Historical Foundations, Euler, Bernoulli & Parity Obstruction

| # | Topic / Query Focus | Key Findings, Authors & Formulas |
|---:|---|---|
| 1 | Euler (1772, 1775) log-sine & $\zeta(2,1)$ | Euler proved $\int_0^{\pi/2} x\ln(\sin x)\,dx = \frac{7}{16}\zeta(3) - \frac{\pi^2}{8}\ln 2$ (1772) and $\zeta(2,1)=\zeta(3)$ in *Meditationes circa singulare serierum genus* (E477, 1775). |
| 2 | Even $\zeta(2k)$ vs odd $\zeta(2k+1)$ cotangent obstruction | $z\cot z = 1 - 2\sum_{k\ge 1}\zeta(2k)(z/\pi)^{2k}$ has strictly even powers because pole pairing $\frac{1}{z-n}+\frac{1}{z+n}=\frac{2z}{z^2-n^2}$ doubles even powers and cancels odd powers identically (formalized in `Zeta3/ParityObstruction.lean`). |
| 3 | Deligne's conjecture (1979) & non-criticality | An integer $m$ is critical for $\mathbb{Q}(0)$ iff neither $\Gamma(s/2)$ nor $\Gamma((1-s)/2)$ has a pole at $s=m$. For $s=3$, dual $1-s=-2$ hits the pole of $\Gamma(-1)$, so $\zeta(3)$ is **non-critical** and governed by Beilinson/Bloch–Kato rather than Deligne periods. |
| 4 | $\zeta(3) = -4\pi^2\zeta'(-2)$ & Bendersky–Adamchik constants | Functional equation maps trivial zero $\zeta(-2)=0$ to $\zeta'(-2)=-\frac{\zeta(3)}{4\pi^2}$; generalized Glaisher–Kinkelin / Bendersky–Adamchik constant $A_2 = \exp(\zeta(3)/(4\pi^2))$. |
| 5 | Dirichlet $\beta(3) = \pi^3/32$ vs $\zeta(3)$ | Opposite parity: $\beta(2n+1) = \frac{(-1)^n E_{2n}}{2(2n)!}(\pi/2)^{2n+1}$ ($E_2=-1 \implies \beta(3)=\pi^3/32$) is critical at *odd* integers, while even $\beta(2n)$ (Catalan $G=\beta(2)$) are non-critical. |
| 6 | Hurwitz $\zeta(3, a)$ & tetragamma $\psi^{(2)}(a)$ | $\zeta(3,a) = -\frac{1}{2}\psi^{(2)}(a)$; at rational $a=p/q$, reductions always retain $\zeta(3)$ as an irreducible basis element alongside $\pi^3$ and Clausen values. |
| 7 | Clausen functions $\operatorname{Cl}_3(\theta)$ & Bernoulli splines | $\operatorname{Cl}_3(\theta)=\sum \cos(k\theta)/k^3$ gives $\operatorname{Cl}_3(0)=\zeta(3)$, $\operatorname{Cl}_3(\pi)=-\frac{3}{4}\zeta(3)$, while the sine partner $\sum \sin(k\theta)/k^3 = \frac{\pi^2\theta}{6} - \frac{\pi\theta^2}{4} + \frac{\theta^3}{12}$ evaluates to polynomials in $\pi$. |
| 8 | Barnes $G$-function & multiple gamma (Adamchik, Choi, Srivastava) | Exact logarithmic closed forms exist in terms of multiple gamma / multiple cosine functions $\mathcal{C}_3(1/4), \mathcal{C}_3(1/6)$ and Catalan's constant $G$. |
| 9 | Borwein–Bradley–Crandall (2000) survey | *J. Comput. Appl. Math.* 121: definitive taxonomy of computational algorithms for $\zeta(3)$ (WZ acceleration, binary splitting, polylogarithm evaluation). |
| 10 | Van der Poorten (1979) *"A proof that Euler missed"* | *Math. Intelligencer* 1: historical & technical exposition of Apéry's 1978 Marseille-Luminy announcement, Markov–Apéry series, and 3-term recurrence verified with Henri Cohen. |

---

## Cluster 2 (Searches 11–20): Rigorous Closed-Form Theory, Periods, Motives & Transcendence

| # | Topic / Query Focus | Key Findings, Authors & Theorems |
|---:|---|---|
| 11 | Timothy Chow (1999) *"What is a Closed-Form Number?"* | Defines **EL-numbers** $\mathbb{E}$ (smallest subfield of $\mathbb{C}$ closed under $\exp, \log$). Under Schanuel's conjecture + Period Conjecture, $\zeta(3) \notin \mathbb{E}$. |
| 12 | Kontsevich & Zagier (2001) *"Periods"* | $\zeta(3) = \iiint_{[0,1]^3}\frac{dx\,dy\,dz}{1-xyz} = \sum_{n=0}^\infty \frac{1}{(n+1)^3}$ certifies $\zeta(3)$ as a period. Period Conjecture: all identities among periods come from additivity, change of variables, and Stokes' theorem. |
| 13 | Grothendieck's Period Conjecture | Predicts $\pi, \zeta(3), \zeta(5), \zeta(7), \dots$ are **algebraically independent over $\mathbb{Q}$** — ruling out ANY nonzero polynomial relation $P(\pi, \zeta(3)) = 0$. |
| 14 | Francis Brown (2012) *"Mixed Tate motives over $\mathbb{Z}$"* | *Annals of Math.* 175: proves Hoffman's $\{2,3\}$-basis theorem for motivic MZVs; $\zeta^{\mathfrak{m}}(3)$ is algebraically independent of $\zeta^{\mathfrak{m}}(2)$. |
| 15 | Deligne–Goncharov motivic fundamental group | Motivic Lie algebra $\mathfrak{g}^{\mathfrak{m}}$ of $\pi_1^{\mathrm{mot}}(\mathbb{P}^1\setminus\{0,1,\infty\})$ is **free** on odd-weight generators $\sigma_3, \sigma_5, \sigma_7, \dots$; $\sigma_3$ is the weight-3 generator dual to $\zeta^{\mathfrak{m}}(3)$. |
| 16 | Schanuel's & André's Generalized Period Conjecture | André's generalized period conjecture subsumes Schanuel's conjecture (Bertolin) and predicts full algebraic independence of $\pi, e, \log 2, \zeta(3)$. |
| 17 | Status of $\zeta(3)/\pi^3 \notin \mathbb{Q}$ | **OPEN.** No peer-reviewed proof exists even ruling out $\zeta(3) = q\pi^3$ ($q\in\mathbb{Q}$), though universally conjectured irrational and transcendental. |
| 18 | Why Baker, Mahler, and Siegel–Shidlovskii fail on $\zeta(3)$ | Three exact structural blocks: (1) $\zeta(3)$ comes from **$G$-functions** (polylogarithms), not $E$-functions (blocks Siegel–Shidlovskii); (2) no $z\mapsto z^d$ scaling equation (blocks Mahler); (3) not a known linear form in $\log\overline{\mathbb{Q}}$ (blocks Baker). |
| 19 | Nesterenko (1996) $\{\pi, e^\pi, \Gamma(1/4)\}$ & $\zeta(3)$ | Nesterenko requires a closed algebraic ODE system (Ramanujan's $P, Q, R$ for $E_2, E_4, E_6$); $\zeta(3)$'s Eichler integral breaks homogeneous modularity with an inhomogeneous period polynomial. |
| 20 | Mahler classification of $\zeta(3)$ | Finite irrationality measure ($w_1(\zeta(3)) < 4.51388$) rules out Liouville behavior; universally conjectured to be an **$S$-number**, though exact class ($S$ vs $T$ vs $U$) is unproven. |

---

## Cluster 3 (Searches 21–30): Ramanujan Formulas, Modular Forms, Eichler Integrals & String Theory

| # | Topic / Query Focus | Key Findings, Authors & Formulas |
|---:|---|---|
| 21 | Ramanujan's formula for $\zeta(2n+1)$ (Berndt, Grosswald 1972) | Exact modular transformation for $\alpha\beta=\pi^2$; at $n=1$: $\frac{1}{\alpha}(\frac{1}{2}\zeta(3)+S(\alpha))+\frac{1}{\beta}(\frac{1}{2}\zeta(3)+S(\beta))=\frac{\alpha^2+5\alpha\beta+\beta^2}{180}$ (verified to $3\times 10^{-85}$ in `zeta3_lab.py`). |
| 22 | Weight-4 Eisenstein series $E_4(\tau)$ Eichler integral | $L(s, E_4) = \zeta(s)\zeta(s-3)$; period polynomial $P_{E_4}(X) = \int_0^{i\infty} E_4(\tau)(X-\tau)^2\,d\tau$ directly encodes $\zeta(3)$ and $\pi^4$. |
| 23 | Gun–Murty–Rath (2011) *Bull. LMS* | *"Transcendental values of certain Eichler integrals"*: links Ramanujan/Grosswald polynomials $R_{2n+1}(z)$ to Eichler integrals and transcendence criteria for odd zeta values. |
| 24 | Murty–Smyth–Wang / Saradha / Chowla–Selberg | Analysis of zeros of Ramanujan polynomials (many lie on $|z|=1$), Chowla–Selberg CM periods, and linear independence of $L$-values. |
| 25 | Simon Plouffe (1998, 2006, 2022) Lambert formulas | Exact decompositions via $S_3(r)=\sum_{k\ge 1}\frac{1}{k^3(e^{\pi r k}-1)}$: $\zeta(3) = 28 S_3(1) - 37 S_3(2) + 7 S_3(4)$ and $\pi^3 = 720 S_3(1) - 900 S_3(2) + 180 S_3(4)$ (both verified to $10^{-82}$ in `zeta3_lab.py`!). |
| 26 | Zagier (1991), Haberland (1983), Kohnen–Zagier | Generating functions of modular period polynomials via Jacobi theta functions; links cusp-form periods directly to odd zeta and MZV relations. |
| 27 | Mock & quantum modular forms (Zagier 2010, Zwegers) | Non-holomorphic Eichler integrals complete mock modular forms to harmonic weak Maass forms; radial limits at cusps produce quantum modular forms linked to odd zetas. |
| 28 | Real-analytic Eisenstein series $E(z,s)$ & Kronecker limits | Fourier constant term carries $\zeta(2s-1)/\zeta(2s)$; higher Kronecker limit formulas for real quadratic fields produce $\zeta(3)$ and MZVs. |
| 29 | Modular graph functions (D'Hoker, Green, Vanhove, Brown) | Genus-one Type II superstring four-graviton amplitudes expand at the cusp in **single-valued MZVs**, with $\zeta(3)$ as the leading weight-3 threshold correction. |
| 30 | Dixit, Kesarwani, Maji, Gupta (2018–2026) | Modern multi-parameter generalizations of Ramanujan's odd-zeta formula, higher Herglotz functions, and Voronoï summation (continuing into 2026). |

---

## Cluster 4 (Searches 31–40): Apéry, Beukers, Hata, Rhin–Viola, Zudilin & 2026 Records

| # | Topic / Query Focus | Key Findings, Authors & Bounds |
|---:|---|---|
| 31 | Beukers (1979) triple integral & shifted Legendre polynomials | $I_n = \iiint_{[0,1]^3} \frac{(x(1-x)y(1-y)z(1-z))^n}{(1-(1-xy)z)^{n+1}}dx\,dy\,dz = a_n\zeta(3) - b_n$; $0 < |I_n| \le \zeta(3)(\sqrt{2}-1)^{4n}$, beating $d_n^3 \sim e^{3n}$. |
| 32 | Masayoshi Hata (1990, 2000) | Asymmetric Legendre-type polynomials improved $\mu(\zeta(3))$ from $13.41782$ down to $8.830284$ (1990) and **$7.377956$** (2000, *Acta Arith.*). |
| 33 | Rhin & Viola (2001) *"The group structure for $\zeta(3)$"* | *Acta Arith.* 97: discovered the **1920-element Weyl group $\Phi \cong W(D_5)$** (vs 120-element $S_5$ for $\zeta(2)$) acting on the 8-parameter triple integral, achieving **$\mu(\zeta(3)) < 5.513891$**. |
| 34 | **Sept 2026 Record: Giraudin (`arXiv:2609.26980`)** | **David Niedbala Giraudin (Sept 22, 2026)** broke the 25-year Rhin–Viola record to **$\mu(\zeta(3)) < 5.5138800$** via narrow-cone perturbation ($52000\times$ RV parameters shifted by $\le 23$), with self-contained C verifier `verify_zeta3.c`. |
| 35 | Wadim Zudilin well-poised ${}_7F_6, {}_9F_8$ series | Unified Beukers/Rhin–Viola/Rivoal via very-well-poised hypergeometric series and hidden permutation symmetries. |
| 36 | Sorokin (1993–1998) Hermite–Padé & Nikishin systems | Reconstructed Apéry's $\zeta(3)$ approximations and $\pi$-transcendence via multiple orthogonal polynomials for Nikishin systems. |
| 37 | Nesterenko (1996, 2003) Mellin–Barnes integrals | Unified contour Barnes integrals with unit-cube integrals for linear forms in zeta values and Nesterenko's linear independence criterion. |
| 38 | Ball–Rivoal (2001) & Zudilin (2001) | Infinitely many odd $\zeta(2k+1)$ are irrational ($\dim \ge c\log s$); at least one of $\zeta(5), \zeta(7), \zeta(9), \zeta(11)$ is irrational. |
| 39 | Fischler–Sprang–Zudilin (2019), Lai–Yu (2020), **Fischler (JLMS 2026)** | FSZ broke $\log s$ for *count* of irrational odd zetas; **Fischler (published 2026)** broke the $\log s$ barrier for *linear independence dimension* using **Siegel's lemma + generalized Shidlovsky multiplicity**, proving $\dim_{\mathbb{Q}}\operatorname{Span}(1,\zeta(3),\dots,\zeta(s)) \ge 0.21\sqrt{s/\log s}$. |
| 40 | Marc Prévost (1996) Remainder Padé Approximants (RPA) | Diagonal $[n/n]$ Padé approximants to the divergent asymptotic expansion of $\zeta(3) - \sum_{k=1}^{n-1} k^{-3}$ reproduce Apéry's sequences organically. |

---

## Cluster 5 (Searches 41–50): Calegari–Dimitrov–Tang (CDT), Bost–Charles, $G$-Functions & K3/CY Geometry

| # | Topic / Query Focus | Key Findings, Authors & Breakthroughs |
|---:|---|---|
| 41 | Calegari–Dimitrov–Tang (JAMS 2025, `arXiv:2408.15403`, ICM 2026 `2510.04156`) | Proved Unbounded Denominators Conjecture, $\mathbb{Q}$-linear independence of $\{1, \zeta(2), L(2,\chi_{-3})\}$, and effective irrationality measures for $L(2,\chi_{-3})$, $\zeta_2(5)$, and $\zeta_5(3)$ via **arithmetic holonomy bounds**. |
| 42 | CDT obstruction for real $\zeta(3)$ and $\pi^3$ | CDT's exact diagnosis: overconvergent modular curves provide clean ODEs for **$p$-adic** odd zetas ($\zeta_2(5), \zeta_5(3)$), whereas **real** odd zetas currently lack an analogous global archimedean ODE over $\mathbb{P}^1\setminus S$ with matching conformal capacity. |
| 43 | Bost slope method & François Charles | Translates auxiliary-polynomial transcendence into Arakelov geometry (Hermitian vector bundles, Harder–Narasimhan slopes, arithmetic surfaces). |
| 44 | Dimitrov (2019), André, Bost–Charles algebraicity criteria | Quantitative arithmetic holonomicity: tension between Liouville denominator bounds and Cauchy/conformal analytic size forces algebraicity or dimension bounds. |
| 45 | Siegel $G$-functions (André, Fischler–Rivoal 2014–2026) | Polylogarithms and odd zeta values belong to the $G$-function world, where a full Siegel–Shidlovskii analogue remains open. |
| 46 | Beukers (1985, 1987) modular parametrization of Apéry's recurrence | Apéry's $\zeta(3)$ differential operator is the Picard–Fuchs equation of a family of elliptic/K3 surfaces parametrized by modular forms on $\Gamma_0(6)^+$. |
| 47 | Golyshev, Zagier, Almkvist–Zudilin Calabi–Yau operators | Order-4 Fuchsian Calabi–Yau operators exhibit integral mirror maps, $p$-adic Frobenius structures, and "Apéry limits" converging to $\zeta(3)$. |
| 48 | Beukers & Peters (1984) K3 surface family | Explicitly identified the one-parameter family of K3 surfaces whose holomorphic 2-form periods satisfy Apéry's differential equation for $\zeta(3)$. |
| 49 | Vasilyev (1996), Zudilin (2003), Krattenthaler–Rivoal | Higher-dimensional unit-cube integrals generalizing Beukers' triple integral to $\zeta(5), \zeta(7), \dots$, proved by Zudilin (2003). |
| 50 | Exact verification of Giraudin (`arXiv:2609.26980`, Sept 2026) | Confirmed exact abstract: 8 parameters = $52000 \times$ Rhin–Viola's optimal tuple shifted by $\le 23$, proving Rhin–Viola's grid optimum is locally suboptimal. |

---

## Cluster 6 (Searches 51–60): WZ Acceleration, Ramanujan Machine, Conservative Matrix Fields & PSLQ

| # | Topic / Query Focus | Key Findings, Authors & Formulas |
|---:|---|---|
| 51 | Markov (1890), Hjortnaes (1953), Amdeberhan–Zeilberger (1997) | Markov anticipated WZ acceleration in 1890; Amdeberhan–Zeilberger generated infinite WZ families for $\zeta(3)$, including $\frac{1}{64}\sum (-1)^k(205k^2+250k+77)(k!)^{10}/((2k+1)!)^5$. |
| 52 | Wedeniwski (1998), Sebah–Gourdon, `y-cruncher` | Wedeniwski's ~5.04 digits/term WZ formula (verified exact to $10^{-81}$ in `zeta3_lab.py`) remains the standard engine in `y-cruncher` for trillion-digit $\zeta(3)$ records. |
| 53 | Pilehrood, Guillera, Campbell, Zudilin (2010–2026) | Markov–WZ pairs and Guillera's "upside-down" duals transform divergent Ramanujan $1/\pi^2$ series into convergent Apéry-like series for $\zeta(3)$. |
| 54 | Zhi-Wei Sun (2015–2026) supercongruences & hypergeometric series | Extensive catalogs of $p^3$–$p^5$ supercongruences and inverse-binomial/harmonic series targeting $\zeta(3)$ and $\pi^3$. |
| 55 | Broadhurst (1998) binary BBP polylogarithmic ladders | Base-$2^{12}$ polylogarithmic ladders enabling BBP-style binary digit extraction for $\zeta(3)$ and $\zeta(5)$. |
| 56 | Ramanujan Machine (Raayoni et al., *Nature* 2021) | Automated MITM/gradient-descent discovery of polynomial continued fractions such as $\frac{8}{7\zeta(3)} = 1 - \frac{1^6}{3\cdot 7 - \frac{2^6}{5\cdot 19 - \dots}}$ (verified to $10^{-84}$ in `zeta3_lab.py`!). |
| 57 | **Conservative Matrix Fields (CMFs)** (Elimelech–Kaminer *PNAS* 2024; Wein et al. `arXiv:2507.08138`, 2025) | Discrete flat connections $M_v(x+w)M_w(x)=M_w(x+v)M_v(x)$ unifying Apéry's proof, Ramanujan Machine CFs, and continuous asymptotics for $\zeta(3)$. |
| 58 | Nesterenko (1996) & Bowman–McLaughlin continued fractions | Systematic theory of polynomial continued fractions $K(a_n/b_n)$ and alternative irrationality CFs for $\zeta(3)$. |
| 59 | Cohen–Rodriguez Villegas–Zagier (2000) | Chebyshev-polynomial alternating-series acceleration ($\sim (3+\sqrt{8})^{-n}$ error), verified to $2.4\times 10^{-83}$ on $\zeta(3)$ in `zeta3_lab.py`. |
| 60 | PSLQ searches (Bailey, Borwein, Broadhurst) | Strict digit-budget & height-bound discipline; large-scale PSLQ hunts have found zero relations between $\zeta(3)$ and powers of $\pi$ up to huge heights, while discovering hundreds of MZV/Feynman/polylog identities. |

---

## Cluster 7 (Searches 61–70): Multiple Zeta Values, Polylogarithms, QFT & Hyperbolic Geometry

| # | Topic / Query Focus | Key Findings, Authors & Theorems |
|---:|---|---|
| 61 | Zagier dimension conjecture & Double Shuffle (IKZ) | $d_k = d_{k-2} + d_{k-3}$ ($\sum d_k x^k = \frac{1}{1-x^2-x^3}$); weight 3 has $d_3 = 1$ because Euler's $\zeta(2,1) = \zeta(3)$ collapses depth 2 to depth 1. |
| 62 | Hoffman's $\{2,3\}$-basis | $\pi^2$ (weight 2) and $\zeta(3)$ (weight 3) are the two atomic generators of the entire MZV space. |
| 63 | Trilogarithm $\operatorname{Li}_3$ functional equations (Lewin, Zagier, Goncharov) | Goncharov's 22-term functional equation for single-valued $\mathcal{L}_3(x)$ governs all weight-3 trilogarithm identities (analogous to Abel's 5-term relation for $\operatorname{Li}_2$). |
| 64 | Zagier's Conjecture on $\zeta_F(3)$ — Goncharov (1995) | Proved for ALL number fields at weight 3 via Goncharov's motivic complex $\Gamma_F(3)$ and higher Bloch group $\mathcal{B}_3(F)$. |
| 65 | QED electron $g-2$ & Feynman periods (Broadhurst, Kreimer, Schnetz, Brown) | $\zeta(3)$ first enters electron $g-2$ at 3-loop QED ($C_3$); Brown proved massless Feynman periods are mixed Tate periods (MZVs) until elliptic/K3/CY thresholds appear at higher loops (Schnetz). |
| 66 | Cosmic Galois group & single-valued MZVs (Brown 2014) | Under the single-valued map, **all even $\zeta^{\mathrm{sv}}(2k) = 0$** while **$\zeta^{\mathrm{sv}}(3) = 2\zeta(3) \neq 0$** survives — cleanly separating $\zeta(3)$ from powers of $\pi$! |
| 67 | Alternating Euler sums & colored MZVs (levels 2, 4, 6) | Level-2 CMZVs yield Landen's $\operatorname{Li}_3(1/2) = \frac{7}{8}\zeta(3) - \frac{\pi^2}{12}\ln 2 + \frac{1}{6}(\ln 2)^3$ (verified to $10^{-82}$). |
| 68 | Multiple $t$-values (Hoffman) & $T$-values (Kaneko–Tsumura) | Odd-index / level-2 variants ($t(3) = \frac{7}{8}\zeta(3)$) with stuffle ($t$) and shuffle ($T$) algebras yielding clean evaluations into $\pi$ and $\zeta(3)$. |
| 69 | Kam Cheong Au (2020–2025) & Charlton (2024) | Systematic WZ-seeds + iterated integrals + `MultipleZetaValues` package proving dozens of Sun/Guillera/Zhao conjectures for $\zeta(3), \zeta(5), \zeta(5,3)$. |
| 70 | Hyperbolic 3-manifolds, Bloch group $\mathcal{B}_2$, Neumann–Zagier | Volumes live in $\mathcal{B}_2$ (dilogarithm); trilogarithm $\mathcal{B}_3$ lifts to higher regulators and $\zeta_F(3)$. |

---

## Cluster 8 (Searches 71–80): Algebraic $K$-Theory, Bloch–Kato, $p$-Adic $\zeta_p(3)$ & Mahler Measures

| # | Topic / Query Focus | Key Findings, Authors & Theorems |
|---:|---|---|
| 71 | Borel's theorem on $K_{2n-1}(\mathbb{Z})$ | For $\mathbb{Z}$ ($r_1=1, r_2=0$), $\operatorname{rank} K_{4k+1}(\mathbb{Z}) = 1$ ($k\ge 1$) and $K_{4k-1}(\mathbb{Z})$ is finite; thus **$K_5(\mathbb{Z}) \otimes \mathbb{Q} \cong \mathbb{Q}$** with Borel regulator proportional to $\zeta(3)$! |
| 72 | Beilinson, Bloch–Kato & ETNC (Huber–Kings, Burns–Flach) | Connects special values of Tate motives to Selmer group orders via Iwasawa Main Conjecture and cyclotomic Euler systems. |
| 73 | $p$-adic $\zeta_p(3)$ (Kubota–Leopoldt vs Deligne–Soulé–Furusho) | Classical Kubota–Leopoldt $p$-adic zeta interpolates negative integers; positive odd $p$-adic $\zeta_p(3)$ lives via **Coleman integration** and Furusho's $p$-adic MZVs. |
| 74 | Irrationality of $p$-adic $\zeta_p(3)$ (Calegari 2005, Beukers 2008, Sprang/Lai/Lupu 2024–2026) | $\zeta_2(3)$ and $\zeta_3(3)$ are proved irrational; recent 2024–2026 results prove $\gg \sqrt{s/\log s}$ odd $p$-adic zeta values are irrational for every prime $p$. |
| 75 | Smyth (1981) Mahler measure | Exact geometric closed form: $m(1 + x + y + z) = \frac{7}{2\pi^2}\zeta(3)$, i.e., $\zeta(3) = \frac{2\pi^2}{7} m(1 + x + y + z)$. |
| 76 | Brunault–Zudilin (2020) *Many Variations of Mahler Measures* | Three-variable Laurent polynomials defining singular K3 surfaces yield Mahler measures proportional to $L(K3, 3)$ and $\zeta(3)$. |
| 77 | Soulé cyclotomic elements & $H^1_f(\mathbb{Q}, \mathbb{Q}_p(3))$ | 1-dimensional Bloch–Kato Selmer group $H^1_f(\mathbb{Q}, \mathbb{Q}_p(3)) \cong \mathbb{Q}_p$; Bloch–Kato exponential map sends $\zeta_p(3)$ to the canonical crystalline extension. |
| 78 | $q$-analogues $\zeta_q(3)$ (Krattenthaler–Rivoal–Zudilin) | Very-well-poised basic hypergeometric series yield irrationality of $\zeta_q(3)$ at integer reciprocals and real quadratic units (some now Lean 4 formalized!). |
| 79 | Finite MZVs & Kaneko–Zagier conjecture | $\mathcal{A} \cong \mathcal{Z}/(\pi^2\mathcal{Z})$; depth-1 $\zeta_{\mathcal{A}}(3) = 0$ since $\sum_{m=1}^{p-1} m^{-3} \equiv 0 \pmod p$ ($p>3$). |
| 80 | $\operatorname{Ext}^1_{\mathrm{MTM}(\mathbb{Z})}(\mathbb{Q}(0), \mathbb{Q}(3))$ | 1-dimensional over $\mathbb{Q}$, generated by the unique non-split Kummer-type motivic extension whose period matrix has $\zeta(3)$ off-diagonal. |

---

## Cluster 9 (Searches 81–90): Function Field Analogues ($\mathbb{F}_q[T]$) & RH Connections

| # | Topic / Query Focus | Key Findings, Authors & Theorems |
|---:|---|---|
| 81 | Carlitz zeta values (Jing Yu 1991; Chang–Yu 2007) | Over $\mathbb{F}_q[\theta]$, Carlitz zeta values at multiples of $q-1$ are algebraic multiples of $\tilde{\pi}^n$ ("even"), while **all remaining ("odd") $\zeta_C(n)$ are PROVED algebraically independent of $\tilde{\pi}$ and each other!** |
| 82 | ABP Criterion (2004) & why it fails over $\mathbb{Q}$ | Positive characteristic has algebraic Frobenius $\tau(x)=x^q$ + rigid analytic trivialization; $\mathbb{Q}$ lacks Frobenius, blocking direct translation of ABP to $\zeta(3)$. |
| 83 | Papanikolas (2008) Tannakian $t$-motivic Galois theorem | Unconditionally proves $\operatorname{tr.deg} = \dim G_X$ over $\mathbb{F}_q(\theta)$ — the exact proved blueprint of Grothendieck's Period Conjecture! |
| 84 | Does RH imply anything about $\zeta(3)$? | RH itself does **not** yield a closed form or improve $\mu(\zeta(3))$ asymptotically (only improves PNT's sub-leading error in $d_n = \operatorname{lcm}(1..n)$). |
| 85 | Hadamard product at $s = 3$ | $\zeta(3) = \frac{2 e^{3B}}{\sqrt{\pi}} \prod_\rho (1 - 3/\rho) e^{3/\rho}$ connects $\zeta(3)$ globally to nontrivial zeros $\rho = \frac{1}{2}+i\gamma$. |
| 86 | Keiper–Li coefficients $\lambda_n$ | Exact formula $\lambda_3 = 1 + \frac{3}{2}\gamma - 3\gamma^2 + \gamma^3 + 3\gamma\gamma_1 + \frac{3}{2}\gamma_2 - \frac{7}{8}\zeta(3)$ embeds $\zeta(3)$ inside Li's RH equivalence ($\lambda_n > 0 \;\forall n$). |
| 87 | Alpöge–Claude (Nov 2026) Lean 4 RH proof | Verified by Lean FRO/Mathlib/con-leche/con-ron/nanoda; uses Stepanov-over-$\mathbb{Z}$ via chains of $\sim\sqrt{X}$ theta lifts on definite unitary groups, Bloch–Kato crystalline extensions, Hörmander $L^2$ $\bar\partial$, and $p$-adic denominator bounds. |
| 88 | Transferability of RH proof machinery to $\zeta(3)$ | Direct transfer to $\zeta(3)$ is open, but 5 specific modules inside `riemann-hypothesis` (`L2DbarEstimates`, `GeometryOfNumbers`, `BlochKatoExtensions`, `EisensteinCokernelLengths`, `ThetaQuasiPeriods`) match modern CDT / Bost–Charles / Rhin–Viola transcendence tools (see `docs/RH_PROOF_CLUES.md`). |
| 89 | RH vs $\operatorname{lcm}(1..n)$ in Beukers/Rhin–Viola | $\ln d_n = \psi(n) = n + O(\sqrt{n}\log^2 n)$ under RH vs $n(1+o(1))$ unconditionally — leading coefficient $c_1=3$ is unchanged, so asymptotic $\mu(\zeta(3))$ is unchanged by RH alone. |
| 90 | Automorphic periods & Rallis inner product formula | Links Petersson norms of theta lifts directly to special $L$-values at critical/near-central points (Shimura, Garrett, Harris–Kudla). |

---

## Cluster 10 (Searches 91–100): Formal Verification (Lean 4/Coq/Isabelle), AI Discovery & Open-Problem Surveys

| # | Topic / Query Focus | Key Findings, Authors & Tools |
|---:|---|---|
| 91 | Liu–Zhang–Zhi (2025, `arXiv:2503.07625`) Lean 4 proof of $\zeta(3) \notin \mathbb{Q}$ | First full Lean 4 formalization (`ahhwuhu/zeta_3_irrational`) via Beukers' triple integral, shifted Legendre polynomials, and *PrimeNumberTheoremAnd*. |
| 92 | Chyzak–Mahboubi–Sibut-Pinote–Tassi (2014) Coq formalization | Certified Apéry's original 1978 proof in Coq using a posteriori verification of creative-telescoping certificates generated in Maple. |
| 93 | Manuel Eberl (2019–2020) Isabelle/HOL formalization | AFP entry `Zeta_3_Irrational` formalizing Beukers' proof via Filaseta's notes + AFP's Prime Number Theorem. |
| 94 | Mathlib4 `riemannZeta` status (2025–2026) | `Mathlib.NumberTheory.LSeries.RiemannZeta` & `HurwitzZetaValues` cover $\zeta(2k)$ via Bernoulli numbers and `RiemannHypothesis`; odd zeta arithmetic remains active research territory. |
| 95 | AI-assisted math discovery (2024–2026) | AlphaProof, AlphaEvolve, FunSearch, Ramanujan Machine, and Claude 2026 demonstrate the winning loop: **automated hypothesis/parameter search + formal Lean 4 certification**. |
| 96 | Creative telescoping (Zeilberger, Chyzak, Koutschan) | Ore-algebra holonomic telescoping (`HolonomicFunctions`) automates Apéry-like recurrences and WZ proofs. |
| 97 | Computer search for Rhin–Viola / Zudilin integrals | Confirmed Rhin–Viola (2001) held for 25 years until Giraudin's Sept 2026 narrow-cone perturbation (`5.5138800`), showing high-scale shifted cones remain unexplored! |
| 98 | Inverse Symbolic Calculator & OEIS A002117 / A013631 | Exhaustive symbolic lookup confirms no elementary closed form exists for $\zeta(3)$ (`A002117`); CF is `A013631` = `[1, 4, 1, 18, 1, 1, 1, 4, ...]`. |
| 99 | Waldschmidt *Open Diophantine Problems* (2004) | Canonical roadmap of folklore conjecture (algebraic independence of $\pi, \zeta(3), \zeta(5), \dots$), linear-form criteria, and saddle-point/denominator tension. |
| 100 | Waldschmidt surveys on Periods & Transcendence (2020–2026) | Synthesizes Kontsevich–Zagier periods, Chudnovsky/Nesterenko theorems, Schanuel's conjecture, and algebraic-group approaches (Wüstholz, Philippon, Waldschmidt). |

*Bonus Searches 101–102:* Exact verification of Rhin–Viola Theorem 5.1 8-parameter integral constraints ($h+m=k+r$, $j+q=l+s$) and exact `arXiv:2609.26980` abstract ($52000\times$ RV parameters shifted by $\le 23$).
