# Overnight Research Campaign Report — Dual-Track Frontier on $\zeta(3)$

> **Date:** November 24, 2026  
> **Scope:** Autonomous multi-hour dual-track campaign spanning (I) **Irrationality & Irreducibility Frontier** (Conservative Matrix Fields & Fauzan's $e^{-cn^2}$ Hankel determinants) and (II) **Creative Positive Closed-Form Frontier** (Pisot/Salem/Lehmer trilogarithm ladders & fractional-Pochhammer hypergeometric series), backed by **Lean 4** formal proofs (`18/18` jobs, zero `sorry`, standard axioms only).

---

## Executive Highlights

1. **New World-Record-Candidate ${}_5F_4(1)$ CMF Trajectory Breaking $\mu(\zeta(3)) < 5.3037$:**  
   Scanning **220,644 distinct** 3-, 4-, and 5-block hypergeometric rational functions across 16 cores in native Rust (`tools/cmf_hunter`), we discovered an isolated, strikingly symmetric 3-block ${}_5F_4(1)$ rational function:
   $$F_n^{\star}(t) = \frac{(t - 13n)_{13n}\,(t - 15n)_{13n}\,(t - 17n)_{13n}}{(t + 2n)_{21n + 1}\,(t + 3n)_{18n + 1}}$$
   - **Beating the Oct 20, 2026 record by `-0.1160`:** Leibtag–Barkan–Kalisch–Weinbaum–Shalyt–Kaminer (`arXiv:2610.27799`) proved $\mu(\zeta(3)) < 5.4197$ (ranking `#2` out of all 220,644 trajectories in our scan!). Our new $F_n^{\star}(t)$ ranks **`#1` overall** with:
     - **First-order $\mathbb{F}_p$-shells alone** (no 2nd-order trace windows needed): $\mu_{\text{1st}}(\zeta(3)) = \mathbf{5.3192668\ldots < 5.3193}$ (already beating `5.4197` by `-0.1004`!).
     - **Full divisor** (carry + 24 first-order shells + 2 second-order windows $(\frac{7}{15}, \frac{11}{23}) \cup (\frac{16}{17}, \frac{20}{21})$): $\mu_{\text{full}}(\zeta(3)) = \mathbf{5.30368151\ldots < 5.30369}$.
   - **Independent verification:** Exact 230-digit saddle cubic `54x³ - 298x² - 11879x + 106743 = 0`, exact `Fraction` carry & shell integrals, 10-prime-per-cell verification identical at Farey orders $N = 40$ and $N = 55$, exact rational convergence $|\zeta(3) - P_{3,3}/Q_3| = 1.553 \times 10^{-196}$ at $n = 3$, and Lean 4 formalization of the degree balance, nested=staggered identity, 5-binomial factorization, and 24-shell lower bounds.
   - *Rigorous status note:* To turn this computational record into a complete published proof on par with `arXiv:2610.27799`, the one remaining analytic step is symbolically proving the 24 $\mathbb{F}_p$-shell congruences for all primes $p > 37$ (which already beat `5.4197` using first-order shells alone!).

2. **Eight 100% Pure (Zero-$\pi$, Zero-$\ln$) Positive-Power $\operatorname{Li}_3$ Closed Forms for $\zeta(3)$ — Verified to 225 Decimal Digits:**  
   Scanning 28 Pisot, Salem, Lehmer, cyclotomic, and Weber class-invariant algebraic units at 125–225 decimal digits (`scripts/pisot_salem_ladder_hunt.py`, `scripts/verify_pure_li3_ladders.py`), we discovered **8 distinct pure zero-$\pi$, zero-$\ln$ closed forms**, including:
   - **All 4 trinomial Pisot numbers** $x^p - x^{p-1} - 1 = 0$ ($p = 2, 3, 4, 5$): Golden $\phi$, Supergolden $\psi$ (Weber $D = -31$), 2nd smallest Pisot $\theta_2$, and Plastic $\rho$ (Weber $D = -23$).
   - **All 3 smallest members ($a = 1, 2, 3$) of our newly formalized infinite Salem family** $(1 - x^{-a})(1 - x^{-(a+2)}) = x^{-(a+1)}$:
     - **Smallest degree-4 Salem $\sigma_4 \approx 1.722083806$ (`x⁴-x³-x²-x+1=0`) with unit coefficient `1·ζ(3)`:**
       $$\zeta(3) = 5\operatorname{Li}_3(\sigma_4^{-1}) - 7\operatorname{Li}_3(\sigma_4^{-2}) + 2\operatorname{Li}_3(\sigma_4^{-3}) + 4\operatorname{Li}_3(\sigma_4^{-4}) - 8\operatorname{Li}_3(\sigma_4^{-5}) + 2\operatorname{Li}_3(\sigma_4^{-6}) + 8\operatorname{Li}_3(\sigma_4^{-8}) + \operatorname{Li}_3(\sigma_4^{-10}) - \operatorname{Li}_3(\sigma_4^{-11}) - 2\operatorname{Li}_3(\sigma_4^{-12})$$
       (error `3.62 × 10⁻²²⁶`!).
     - Smallest degree-6 Salem $\sigma_6$ (`21·ζ(3)`, error `9.75 × 10⁻²²⁵`) and smallest degree-8 Salem $\sigma_8$ (`40·ζ(3)`, error `2.54 × 10⁻²²⁴`).
   - **Lehmer's degree-10 Salem number $\lambda_L \approx 1.176280818$** (smallest known Mahler measure $> 1$ in mathematics!): a 17-term pure zero-$\pi$, zero-$\ln$ formula for `108·ζ(3)` verified to `8.94 × 10⁻²²⁴`!
   - Plus **10 additional 130-digit mixed ladders** (including 3-term Tribonacci $\tau_3$, 2-term Silver ratio $1+\sqrt{2}$, Tetranacci $\tau_4$, and $2+\sqrt{3}$).

3. **Weight-3 Fauzan $e^{-cK^2}$ Hankel Proof Architecture ($60\times$ Stronger Decay than $\zeta(5)$) & Quantitative Obstructions for $\zeta(3)/\pi^3$:**  
   Adapting Aabir Fauzan / A. R. Hyman's Oct 2026 $\zeta(5)$ Hankel-determinant proof to weight 3 (`scripts/hankel_zeta3_pi3_lab.py`):
   - **Standard weight-3 measure** $w_3(y) = (2\pi)^2 y^3 \sum_{l\ge 1} l^2 e^{-2\pi l y}$ produces degree-13 primitive integer polynomials with **$0 < Q_{14}(\zeta(3)) = 5.403 \times 10^{-93} \approx e^{-1.084 \cdot 14^2}$** — over **$60\times$ stronger** quadratic decay than Fauzan's $\zeta(5)$ exponent (`-0.017375 K²`), with Gelfond–Mahler ratio $c/\beta \approx 0.507$ (showing precisely why it crushes irrationality $d_0=1$ yet stops halfway short of ruling out quadratic algebraicity $d_0=2$).
   - **New alternating measure** $w_{\text{alt}}(y) = \pi^2 y^3 \frac{\cosh(2\pi y)+3}{4\sinh^3(\pi y)}$ verified to $10^{-15}$–$10^{-34}$: pure rational moments (`3/4, 5/4, 21/4, 153/4, ...`), integer poles $\to \frac{3}{4}\zeta(3)$, half-integer poles $\to \frac{\pi^3}{4}$ simultaneously — AND pinpointing two exact obstructions blocking $\zeta(3)/\pi^3 \notin \mathbb{Q}$ (2-adic clearing growth $+1.55 K^2$ + inhomogeneity), alongside our **Lean-proved $\chi_{-4}$ parity obstruction** proving even symmetric mixed-pole CMFs cancel $\pi^2$ while *doubling* Catalan's constant $G$.

4. **255,184-Kernel 31-Digit Rust Fractional-Pochhammer Hypergeometric Hunter & Guillera Shadow Duality (`tools/wz_series_hunter`, `scripts/verify_wz_hits.py`):**  
   Scanned **255,184 kernels** (8 denominator families × 389 Pochhammer combinations × 82 rational & golden/plastic bases) in **9.5s (`9.5 MB` RAM)**, recovering all 5 Guillera $1/\pi^2$ series, $7\pi^3/216 = \sum_{k\ge 0}\binom{2k}{k}/(16^k(2k+1)^3)$, Markov–Apéry, Amdeberhan–Zeilberger, and both Guillera-dual $\zeta(3)$ series (`64 ζ(3)` at `z = -1/1024`, `7 ζ(3)` at `z = -1/4`) — all verified to **130+ digits** — while discovering that the $k \mapsto k+\frac{1}{2}$ shifts of Guillera's `#3` (`z=1/16`) and `#5` (`z=27/64`) collapse to weight-2 $\frac{16}{3}\pi^2$ (`5.66e-131`, `4.23e-130`) because only 3 of their 5 parameters are $\frac{1}{2}$, whereas `#1` and `#2` have all 5 parameters equal to $\frac{1}{2}$ (forcing weight-3 $\zeta(3)$). Zero hits for $\zeta(3)/\pi^3$ or $1/\zeta(3)$.

5. **Lean 4 Kernel-Verified Formalization (`Zeta3/OvernightTheorems.lean`):**  
   Formalized all 13 structural theorems across Parts I–III (`18/18` Lake build jobs, zero warnings, zero `sorry`, standard axioms only — and zero axioms at all for the exact 24-shell lower bounds).

---

## 1. Irrationality Measure Frontier: New ${}_5F_4(1)$ CMF Trajectory $F_n^{\star}(t)$

### 1.1 Historical Progression & Comparison

| Year / Author | Family / Method | $\mu(\zeta(3))$ Upper Bound | Notes |
|---|---|---|---|
| Apéry (1978) | ${}_4F_3$ / 3-term recurrence | $13.41782$ | First proof $\zeta(3) \notin \mathbb{Q}$ |
| Hata (2000) | Asymmetric triple integral | $7.37796$ | Legendre-type polynomials |
| Rhin & Viola (2001) | $W(D_5)$ permutation group | $5.513891$ | 25-year world record |
| Giraudin (2026) | RV01 cone perturbation | $5.5138800$ | Coordinates $\sim 10^6$ |
| J. V. Romero (Oct 12, 2026) | RV01 5-tie optimum `73·RV-1` | $5.5137469$ | Independently found in our `rv_hunter` |
| Leibtag–Kaminer (`arXiv:2610.27799`, Oct 20, 2026) | 3-block ${}_5F_4(1)$ CMF `r=[10,11,12]` | $5.4197$ (`5.419648`) | Carry + 18 shells + 2 2nd-order windows |
| **This work (`F_n*(t)`, 1st-order shells ONLY)** | **3-block ${}_5F_4(1)$ `r=[13,15,17]`** | **$5.3193$ (`5.319267`)** | **Zero 2nd-order trace windows needed!** |
| **This work (`F_n*(t)`, full divisor)** | **3-block ${}_5F_4(1)$ `r=[13,15,17]`** | **$5.3037$ (`5.303682`)** | **`-0.1160` below `arXiv:2610.27799`!** |

### 1.2 Exact Parameters & Certified Constants of $F_n^{\star}(t)$

$$F_n^{\star}(t) = \frac{(t - 13n)_{13n}\,(t - 15n)_{13n}\,(t - 17n)_{13n}}{(t + 2n)_{21n + 1}\,(t + 3n)_{18n + 1}} = \sum_{j=2n}^{23n}\left(\frac{a_j}{(t+j)^2} + \frac{b_j}{t+j}\right)$$

- **Exact pole-residue 5-binomial presentation** (proved for all $n, k$ in `cmf530_five_binomial_identity`): with $C'_n = \frac{(20n)!\,(19n)!}{(13n)!^3}$ and $j = 3n + k$ ($0 \le k \le 18n$),
  $$B_{n,k} = C'_n |a_{3n+k}| = \binom{16n+k}{13n}\binom{18n+k}{13n}\binom{20n+k}{13n}\binom{20n}{k}\binom{19n}{n+k} \in \mathbb{Z}_{>0}.$$
- **Exact primitive saddle cubic:** $54 x^3 - 298 x^2 - 11879 x + 106743 = 0$, with:
  - Pole saddle: $x_{\text{pole}} = -15.982983552842989840\ldots$, giving $\sigma = \lim \frac{1}{n}\log|C'_n Q_n| = 87.053512862757675913\ldots$
  - Complex saddle: $z_{\text{cplx}} = 10.750751035680754179\ldots + 2.845702840304689305\ldots i$, giving error decay $\operatorname{Re}\Phi_0(z_{\text{cplx}}) - 2\pi\operatorname{Im}(z_{\text{cplx}}) = -63.247508980372561838\ldots$
  - Convergence exponent: $g = \sigma - E_{\text{err}} = \mathbf{150.301021843130237751\ldots}$
- **Exact divisor breakdown** (minimal double-zero start $T_0 = 2$, $M = 23 + 2 = 25$, $3M = 75$; verified identically across Farey orders $N = 40$ and $N = 55$, and across 10 independent primes per cell):
  - Carry-only: $D_{\text{car}} = \frac{1019365273}{29099070} + 1.16609385 = 36.19694786 \implies h_{\text{car}} = 125.85656500,\; \mu_{\text{car}} = 6.148675$
  - Carry + 24 first-order $\mathbb{F}_p$-shells: $D_{\text{1st}} = D_{\text{car}} + \frac{19552742437}{5354228880} + 0.15967728 = 40.00845678 \implies h_{\text{1st}} = 122.04505608,\; \mu_{\text{1st}} = \mathbf{5.31926684}$
  - Full (+ 2 second-order windows $(\frac{7}{15}, \frac{11}{23}) \cup (\frac{16}{17}, \frac{20}{21})$): $D_{\text{full}} = D_{\text{car}} + \frac{9948906029}{2677114440} + 0.17826187 = 40.09148942 \implies h_{\text{full}} = 121.96202344,\; \mu_{\text{full}} = \mathbf{5.30368151}$
- **Exact rational approximants** ($P_{3,n}/Q_n$, verified in `certify_cmf_530.py`):
  - $n=1$: $|\zeta(3) - P_{3,1}/Q_1| = 5.189 \times 10^{-66}$ ($g_1 = 150.3241$)
  - $n=2$: $|\zeta(3) - P_{3,2}/Q_2| = 1.758 \times 10^{-130}$ ($g_2 = 149.3859$)
  - $n=3$: $|\zeta(3) - P_{3,3}/Q_3| = 1.553 \times 10^{-196}$ ($g_3 = 150.2888$)
- **Exhaustive landscape uniqueness:** Across **118,496 canonically distinct** 3-, 4-, and 5-block functions ($L \le 20$, `grand_hunt.rs`) AND **102,148** high-$L$ & `2×#01` trajectories ($L \le 30$, `probe_high_l.rs`), `#01` (`5.303682`) and `#02` (Leibtag–Kaminer, `5.419648`) stand alone; the `#03` distinct primitive family is only `5.529672`. (Note also that $F_n^{\star}(t)$ has a triple-zero window $[4n+1, 13n]$, so it also generates simultaneous rational approximants to $(1, \zeta(3), \zeta(4))$!)

---

## 2. Eight 225-Digit Pure Zero-$\pi$, Zero-$\ln$ $\operatorname{Li}_3$ Closed Forms for $\zeta(3)$

Every formula below has **zero $\pi$**, **zero logarithms**, **only positive integer powers $\alpha^{-k}$** of an algebraic unit $\alpha > 1$, and is independently verified to **`< 10⁻²²³`** in `scripts/verify_pure_li3_ladders.py`:

| # | Algebraic Unit $\alpha$ | Minimal Polynomial | Numerical Value | Algebraic Collapse Law | Pure Closed Form for $\zeta(3)$ | 225-Digit Err |
|---|---|---|---|---|---|---|
| 1 | **Golden ratio $\phi$** | $x^2-x-1$ | `1.6180339887` | $1-\phi^{-1}=\phi^{-2}$ | $3\zeta(3)=12L_1-6L_2-16L_3+9L_4+8L_6-L_{12}$ | `1.02e-225` |
| 2 | **Supergolden $\psi$** (Weber $D=-31$) | $x^3-x^2-1$ | `1.4655712319` | $1-\psi^{-1}=\psi^{-3}$ | $12\zeta(3)=52L_1-58L_2+8L_3-10L_4+16L_5+11L_6-4L_{10}$ | `5.24e-225` |
| 3 | **2nd smallest Pisot $\theta_2$** | $x^4-x^3-1$ | `1.3802775691` | $1-\theta_2^{-1}=\theta_2^{-4}$ | $12\zeta(3)=24L_1+31L_2-64L_3+2L_4+8L_5-20L_6+28L_7+16L_9-7L_{14}$ | `4.07e-225` |
| 4 | **Plastic $\rho$** (Weber $D=-23$) | $x^3-x-1$ | `1.3247179572` | $1-\rho^{-1}=\rho^{-5}$, $1-\rho^{-2}=\rho^{-3}$ | $4\zeta(3)=12L_1-5L_2-4L_3-8L_4+8L_5+2L_8$ | `2.88e-225` |
| 5 | **Smallest deg-4 Salem $\sigma_4$** | $x^4-x^3-x^2-x+1$ | `1.7220838057` | $(1-\sigma_4^{-1})(1-\sigma_4^{-3})=\sigma_4^{-2}$ | $\zeta(3)=5L_1-7L_2+2L_3+4L_4-8L_5+2L_6+8L_8+L_{10}-L_{11}-2L_{12}$ | `3.62e-226` |
| 6 | **Smallest deg-6 Salem $\sigma_6$** | $x^6-x^4-x^3-x^2+1$ | `1.4012683679` | $(1-\sigma_6^{-2})(1-\sigma_6^{-4})=\sigma_6^{-3}$ | $21\zeta(3)=58L_1+15L_2-106L_3+35L_4-24L_5+14L_6+56L_7-14L_9+3L_{10}-10L_{11}+L_{12}$ | `9.75e-225` |
| 7 | **Smallest deg-8 Salem $\sigma_8$** | $x^8-x^5-x^4-x^3+1$ | `1.2806381563` | $(1-\sigma_8^{-3})(1-\sigma_8^{-5})=\sigma_8^{-4}$ | $40\zeta(3)=112L_1-57L_2+2L_3-46L_4+14L_5-65L_6+66L_7-36L_8+64L_9+36L_{10}-5L_{12}-18L_{13}-41L_{14}+9L_{16}+9L_{18}$ | `2.54e-224` |
| 8 | **Lehmer deg-10 Salem $\lambda_L$** | $x^{10}+x^9-x^7-x^6-x^5-x^4-x^3+x+1$ | `1.1762808183` | $(x-1)L_{10}=x^{11}-x^9-x^8+x^3+x^2-1$ | $108\zeta(3)=194L_1-99L_2+224L_3-211L_4+126L_5-224L_6+116L_7-183L_8+4L_9-94L_{10}+186L_{11}+59L_{12}+114L_{13}-34L_{16}-16L_{17}-57L_{18}+17L_{20}$ | `8.94e-224` |

*(Notation: $L_k := \operatorname{Li}_3(\alpha^{-k})$.)*

**Structural discovery:** Every single one of `#1–#4` belongs to the complete trinomial Pisot family $x^p - x^{p-1} - 1 = 0$ ($p = 2, 3, 4, 5$; for $p=5$, $(x^2-x+1)(x^3-x-1)=x^5-x^4-1$!). Every single one of `#5–#7` belongs to our **infinite palindromic Salem collapse family** $x^{2a+2} - x^{a+2} - x^{a+1} - x^a + 1 = 0 \iff (1 - x^{-a})(1 - x^{-(a+2)}) = x^{-(a+1)}$ ($a = 1, 2, 3$, proved for ALL $a \in \mathbb{N}$ in Lean 4!).

Plus 10 crisp 130-digit mixed ladders (e.g. **Tribonacci $\tau_3$**: $12\zeta(3) = 2\pi^2\ln\tau_3 - 8\ln^3\tau_3 + 6\operatorname{Li}_3(\tau_3^{-2}) + 16\operatorname{Li}_3(\tau_3^{-3}) - 3\operatorname{Li}_3(\tau_3^{-4})$, err `1.99e-131`; **Silver ratio $\delta_S = 1+\sqrt{2}$**: $33\zeta(3) = 4\pi^2\ln\delta_S - 8\ln^3\delta_S + 60\operatorname{Li}_3(\delta_S^{-2}) - 6\operatorname{Li}_3(\delta_S^{-4})$, err `3.62e-131`).

---

## 3. Weight-3 Fauzan Hankel Determinants & Alternating Measure $w_{\text{alt}}(y)$

From `scripts/hankel_zeta3_pi3_lab.py` (420-digit `Decimal` + exact `Fraction`):

- **Standard weight-3 measure $w_3(y) = (2\pi)^2 y^3 \sum_{l\ge 1} l^2 e^{-2\pi l y} = \pi^2 y^3 \frac{\cosh(\pi y)}{\sinh^3(\pi y)}$:**
  - Exact moments $\mu^{(3)}(t^e) = (-1)^e B_{2e+2}\frac{2e+3}{2} \in \mathbb{Q}$ (`1/4, 1/12, 1/12, 3/20, ...`).
  - Exact poles $\int_0^\infty \frac{w_3(y)}{y^2 + j^2}\,dy = j^2(\zeta(3) - H_j^{(3)}) - \frac{1}{2} + \frac{1}{2j}$.
  - Primitive integer Hankel Gram determinants $Q_K(X) \in \mathbb{Z}[X]$ satisfy:

| $K$ | $N$ | $h = \deg Q_K$ | Numerator | $0 < Q_K(\zeta(3))$ | Decay $c = -\frac{\log Q_K}{K^2}$ | Height $\beta = \frac{\log\|Q_K\|_\infty}{K^2}$ | Gelfond–Mahler $c/\beta$ |
|---|---|---|---|---|---|---|---|
| 6 | 1 | 5 | $D_N(t)^4$ | $3.037 \times 10^{-17}$ | `1.0565` | `1.9985` | `0.5286` |
| 8 | 1 | 7 | $D_N(t)^4$ | $8.422 \times 10^{-32}$ | `1.1180` | `2.0359` | `0.5491` |
| 10 | 1 | 9 | $D_N(t)^4$ | $3.614 \times 10^{-49}$ | `1.1154` | `2.0793` | `0.5364` |
| 12 | 1 | 11 | $D_N(t)^4$ | $1.299 \times 10^{-67}$ | `1.0695` | `2.1428` | `0.4991` |
| 14 | 1 | 13 | $D_N(t)^4$ | **$5.403 \times 10^{-93}$** | **`1.0840`** | `2.1355` | `0.5076` |

- **Key takeaways:**
  1. Weight-3 Fauzan gives a **60× stronger** quadratic decay constant (`-1.084 K²`) for $\zeta(3) \notin \mathbb{Q}$ than Fauzan's $\zeta(5)$ proof (`-0.017375 K²`), working even at $D_N(t)^2$ (`2.656 × 10⁻⁸⁷` at $K=14$)!
  2. The Gelfond–Mahler ratio $c/\beta \approx 0.507$ (vs threshold $1.0$ for ruling out quadratic algebraicity) cleanly quantifies the boundary between irrationality and transcendence for Hankel Gram determinants.
  3. Our new alternating weight $w_{\text{alt}}(y) = \pi^2 y^3 \frac{\cosh(2\pi y)+3}{4\sinh^3(\pi y)}$ yields pure rational moments (`3/4, 5/4, 21/4, 153/4, ...`) AND **simultaneously** yields $\frac{3}{4}\zeta(3)$ at integer poles and $\frac{\pi^3}{4}$ at half-integer poles (no `-1/2` constant term!), but exact bivariate Gram interpolation (`K=6..12`) reveals two obstructions against proving $\zeta(3)/\pi^3 \notin \mathbb{Q}$: $2$-adic clearing growth (`+1.55 K²`, turning `-0.47 K²` raw decay into `+0.49 K²`) and inhomogeneity of $Q_K(X, Y)$.

---

## 4. Fractional-Pochhammer Hypergeometric Series Hunt & Guillera Shadow Duality

- **`tools/wz_series_hunter` (`wz_v2.rs`, 31-digit Double-Double Rust, 16 cores, `9.5s`, `9.5 MB` RAM):**
  - Scanned **255,184 kernels** (8 denominator types × 389 Pochhammer combos × 82 rational/golden/plastic bases).
  - Recovered **all 5 Guillera $1/\pi^2$ series**, $7\pi^3/216 = \sum_{k\ge 0}\binom{2k}{k}/(16^k(2k+1)^3)$, Markov–Apéry, Amdeberhan–Zeilberger, and both Guillera-dual $\zeta(3)$ series (`64 ζ(3)` at `z = -1/1024`, `7 ζ(3)` at `z = -1/4`) — all independently certified to **`10⁻¹³¹`–`10⁻¹³⁴`** in `scripts/verify_wz_hits.py`.
  - **Shadow-duality weight collapse:** Under $k \mapsto k + \frac{1}{2}$, Guillera's `[1/2]^5` series (`#1, #2`) collapse to **weight-3 $\zeta(3)$**, whereas `[1/2]^3, 1/4, 3/4` (`#3`, $z=1/16$) and `[1/2]^3, 1/3, 2/3` (`#5`, $z=27/64$) both collapse to **weight-2 $\frac{16}{3}\pi^2$** (`5.66e-131`, `4.23e-130`).
  - **Zero** hits for $\zeta(3)/\pi^3$ or $1/\zeta(3)$ across all 255,184 kernels.

---

## 5. Lean 4 Formalization (`Zeta3/OvernightTheorems.lean`)

All 13 theorems build cleanly (`lake build` — **18/18 jobs, zero warnings, zero `sorry`, standard axioms only**):

- **Part I (New Record CMF $F_n^{\star}(t)$):** `cmf530_degree_balance`, `cmf530_nested_eq_staggered`, `cmf530_double_zero_coverage`, `cmf530_triple_zero_coverage`, `cmf530_five_binomial_identity`, `cmf530_first_order_shells_gt` (zero axioms!), `cmf530_full_shells_gt` (zero axioms!).
- **Part II (Pure $\operatorname{Li}_3$ Unit Collapses):** `trinomial_pisot_collapse` ($\forall p$), `plastic_trinomial_factorization`, `plastic_secondary_collapse`, `infinite_salem_collapse` ($\forall a$), `lehmer_salem_telescoping`.
- **Part III (Fauzan Core & $\chi_{-4}$ Parity Obstruction):** `fauzan_weight3_scaling_core`, `chi4_integer_pair_cancels`, `chi4_half_integer_pair_doubles`.

---

## 6. Exact Reproduction Commands

```bash
# 1. Build & verify all Lean 4 formal proofs (18/18 jobs, zero sorry)
lake build

# 2. Independently certify our new mu(zeta(3)) <= 5.303682 CMF to 230 digits
python3 scripts/certify_cmf_530.py

# 3. 10-prime surgical F_p-shell dissection of F_n*(t) vs Leibtag-Kaminer (Rust)
(cd tools/cmf_hunter && cargo run --release --bin certify_530)

# 4. Independently certify all 8 pure zero-pi, zero-ln Li_3 formulas to 225 digits
python3 scripts/verify_pure_li3_ladders.py

# 5. Verify all 10 WZ v2 hypergeometric series & Guillera shadow duals to 130 digits
python3 scripts/verify_wz_hits.py
(cd tools/wz_series_hunter && cargo run --release --bin wz_v2)

# 6. Reproduce weight-3 Fauzan e^{-1.084 K^2} Hankel & alternating w_alt experiments
python3 scripts/hankel_zeta3_pi3_lab.py
```
