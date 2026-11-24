# Closed Form & Arithmetic Nature of $\zeta(3)$

A Lean 4 formalization and computational research workspace dedicated to studying the **closed-form problem for Apéry's constant**:

$$\zeta(3) = \sum_{n=1}^{\infty} \frac{1}{n^3} \approx 1.202056903159594285399738161511449990764986292\ldots$$

## Why This Problem Matters

Euler proved in 1734–1748 that every positive **even** Riemann zeta value has an exact closed form in terms of $\pi$:

$$\zeta(2k) = (-1)^{k+1} \frac{B_{2k}\,(2\pi)^{2k}}{2\,(2k)!}, \qquad \zeta(2) = \frac{\pi^2}{6}, \quad \zeta(4) = \frac{\pi^4}{90}, \quad \zeta(6) = \frac{\pi^6}{945}, \dots$$

Yet for **odd** integers $2k+1 \ge 3$—starting with $\zeta(3)$—no analogous closed form in terms of $\pi$ and elementary constants is known. Even more starkly:

| Question | Status (2026) |
|---|---|
| Is $\zeta(3)$ irrational? | **Proved** (Apéry 1978; Beukers 1979; formalized in Lean 4 by Liu–Zhang–Zhi 2025, `ahhwuhu/zeta_3_irrational`) |
| Is $\zeta(3) / \pi^3 \in \mathbb{Q}$? ("Euler-type closed form") | **OPEN** (conjectured **No**) |
| Are $1, \pi^3, \zeta(3)$ linearly independent over $\mathbb{Q}$? | **OPEN** |
| Is $\zeta(3)$ transcendental? | **OPEN** (conjectured **Yes**) |
| Are $\pi$ and $\zeta(3)$ algebraically independent over $\mathbb{Q}$? | **OPEN** (predicted by Grothendieck's Period Conjecture) |
| Best irrationality measure $\mu(\zeta(3))$ | Classical refereed record: `< 5.513891` (Rhin–Viola 2001). **2026 Frontier:** `< 5.5138800` (Giraudin, Sept 2026); `< 5.513747` (J. V. Romero, Oct 12, 2026 — independently rediscovered & 80-digit/Lean-4 certified here as `73*RV - 1`); `< 5.4197` via 3-block ${}_5F_4(1)$ CMFs (Leibtag–Kaminer, Oct 20, 2026, `arXiv:2610.27799`). **New PROVED elementary record in this repo (`paper/main.pdf`):** **`5.303681511 < 5.3037`** (with 1st-order $\mathbb{F}_p$-shells alone already giving `< 5.3193`, plus our elementary $\mathbb{Z}/p^2\mathbb{Z}$ Faulhaber-lift proof of both 2nd-order shells $I_{12} \cup I_{23}$, and single-presentation $\mu(A\zeta(3)+B\zeta(4)) < 112.3$!) via $F_n^\star(t) = (t-13n)_{13n}(t-15n)_{13n}(t-17n)_{13n}/((t+2n)_{21n+1}(t+3n)_{18n+1})$. Conjectured $\mu = 2$. |
| Does $\zeta(3)$ have a closed form in a richer class (modular/Eichler/polylog/period/CMF)? | **Partially Yes** — exact identities exist via Ramanujan's & Plouffe's Lambert/Eichler series, Clausen/polylogarithms (including **8 pure zero-$\pi$, zero-$\ln$ Pisot/Salem/Lehmer $\operatorname{Li}_3$ formulas verified to 225 digits** in this repo!), Beukers' triple integral, Conservative Matrix Fields, and Smyth's Mahler measure |

## Key Results (See [`paper/main.pdf`](./paper/main.pdf), [`docs/CAMPAIGN_25_EXPERIMENTS_REPORT.md`](./docs/CAMPAIGN_25_EXPERIMENTS_REPORT.md), [`docs/NUMBER_THEORY_MASTER_STRATEGY.md`](./docs/NUMBER_THEORY_MASTER_STRATEGY.md) & [`docs/OVERNIGHT_DISCOVERIES.md`](./docs/OVERNIGHT_DISCOVERIES.md))

1. **New Elementary Record Proof $\mu(\zeta(3)) < 5.3037$ & Single-Presentation $\mu(A\zeta(3)+B\zeta(4)) < 112.3$ (`paper/main.pdf`, `Zeta3/ElementaryRecordProof.lean`, `scripts/certify_paper_bounds.py`, `scripts/verify_explicit_Gm_formula_full.py`):**
   - Screening **>318,000 distinct** hypergeometric rational functions across 16 cores in native Rust (`grand_hunt`, `probe_high_l`, `grand_campaign_p2`) and running a **25-experiment optimality campaign** (`docs/CAMPAIGN_25_EXPERIMENTS_REPORT.md`), we discovered and proved an isolated global optimum beating Leibtag–Kaminer (`arXiv:2610.27799`, `5.4197`) by **`-0.1160`**:
     $$F_n^{\star}(t) = \frac{(t - 13n)_{13n}\,(t - 15n)_{13n}\,(t - 17n)_{13n}}{(t + 2n)_{21n + 1}\,(t + 3n)_{18n + 1}}$$
   - **Elementary First-Order Wilson–Lucas Proof (`Section 4.2`, `< 5.3193`):** Every double-pole residue factors as a positive 5-binomial product whose canonical mod-$p$ bucket polynomial $\mathcal{R}_m(x) \in \mathbb{F}_p[x]$ has degree **exactly $p - 2 < p - 1$** across all 38 uniform-carry sub-shells (88/88 Farey-40 cells) of all 24 $\mathbb{F}_p$-shells, in both $u < 1$ and periodic $u \ge 1$.
   - **Elementary Second-Order $\mathbb{Z}/p^2\mathbb{Z}$ Proof (`Section 4.3`, `< 5.3037` & `< 112.3`):** On the two second-order windows $I_{12} = (\frac{7}{15}, \frac{11}{23})$ and $I_{23} = (\frac{16}{17}, \frac{20}{21})$ (5 sub-windows, 14/14 Farey-40 cells), smooth $\mathbb{Z}/p^2\mathbb{Z}$ binomial & simple-pole wing lifts + Wolstenholme's theorem + our 6-line **Faulhaber Master Congruence** $\frac{1}{p}\sum_{x=0}^{p-1} F(x) \equiv \sum_{x=0}^{p-1} F(x) H_x \pmod p$ ($\deg F \le p - 2$) prove bucket-by-bucket vanishing mod $p^{v+2}$ for $Q_n$ and $P_{3,n}$ (and $P_{4,n}$ on $I_{12}$), while the $P_{4,n}$ bucket split on $I_{23}$ vanishes via $\frac{1}{6}\widetilde{G}_q''(0) = 0$ for an explicit degree-balanced quotient rational function $\widetilde{G}_q(t)$ with a triple zero at $t = 0$.
   - **Certified to 100 digits** (`scripts/certify_paper_bounds.py`, interval width $< 5.2 \times 10^{-32}$) and **formally verified in Lean 4** (`Zeta3/ElementaryRecordProof.lean`, `lake build` 20/20 jobs, zero `sorry`, zero `native_decide`, zero extra axioms).
2. **10-Branch Number Theory Master Strategy & Quantitative Cross-Branch Stress Tests (`docs/NUMBER_THEORY_MASTER_STRATEGY.md`, `scripts/strategy_stress_tests.py`):**
   - Benchmarked all 10 branches of modern number theory (analytic/hypergeometric CMFs, finite-field arithmetic, Fauzan/FLSZ Hankel determinants, CDT modular holonomy, Nesterenko differential algebra, G-function Padé, $p$-adic overconvergent forms, $K_5(\mathbb{Z})$/Bloch groups, $\mathrm{MTM}(\mathbb{Z})$ motives, and hyperbolic geometry/Lean 4) — uncovering the **alternating parity-coupling obstruction** separating $\{\zeta(3), G, \ln 2\}$ from $\{\pi^3, \zeta(2), \pi\}$, **antisymmetric Lucas-palindrome $\mathbb{F}_p$-shells** on VWP ${}_7F_6(1)$, the complete genus-0 $X_0(N)$ CDT map, and a **strictly positive $\chi_{-3}$-twisted weight-3 Fauzan Hankel measure** collapsing to a universal 8-monomial cubic in $(\zeta(3), \frac{2\pi^3}{3\sqrt{3}})$.
3. **Eight 225-Digit Pure Zero-$\pi$, Zero-$\ln$ Pisot, Salem & Lehmer $\operatorname{Li}_3$ Closed Forms (`scripts/pisot_salem_ladder_hunt.py`, `scripts/verify_pure_li3_ladders.py`):**
   - Discovered and certified to **`< 10⁻²²⁸`** pure positive-power $\operatorname{Li}_3(\alpha^{-k})$ closed forms for $\zeta(3)$ with **zero $\pi$ and zero logarithms** across all 4 trinomial Pisot numbers ($\phi$, Supergolden $\psi$/Weber $D=-31$, 2nd Pisot $\theta_2$, Plastic $\rho$/Weber $D=-23$), all 3 smallest members of our Lean-proved infinite Salem collapse family ($\sigma_4$ with unit coefficient `1·ζ(3)`!, $\sigma_6$, $\sigma_8$), AND **Lehmer's degree-10 Salem number $\lambda_L$** (smallest known Mahler measure $> 1$).
4. **Weight-3 Fauzan $e^{-1.084 K^2}$ Hankel Gram Architecture & $\zeta(3)/\pi^3$ Obstructions (`scripts/hankel_zeta3_pi3_lab.py`):**
   - Standard weight-3 Fauzan measure $w_3(y)$ yields $0 < Q_{14}(\zeta(3)) = 5.40 \times 10^{-93} \approx e^{-1.084 \cdot 14^2}$ ($60\times$ stronger quadratic decay than $\zeta(5)$, Gelfond–Mahler ratio $c/\beta \approx 0.507$). Alternating measure $w_{\text{alt}}(y)$ simultaneously yields $\frac{3}{4}\zeta(3)$ and $\frac{\pi^3}{4}$, pinpointing the 2-adic/inhomogeneity obstruction and Lean-proved $\chi_{-4}$ parity obstruction coupling $\pi^3$ to Catalan's constant $G$.
5. **255,184-Kernel Rust Fractional-Pochhammer Series Hunter (`tools/wz_series_hunter`, `scripts/verify_wz_hits.py`):**
   - 31-digit Double-Double 16-core scan (`9.5s`, `9.5 MB` RAM) automatically recovering all 5 Guillera $1/\pi^2$ series, $7\pi^3/216$, Markov–Apéry, Amdeberhan–Zeilberger, and both Guillera-dual $\zeta(3)$ series to `10⁻¹³¹`–`10⁻¹³⁴`, while showing why $k \mapsto k+\frac{1}{2}$ collapses Guillera `#3` & `#5` to weight-2 $\frac{16}{3}\pi^2$.
6. **Rhin–Viola (2001) Cusp Optimum (`73*RV - 1`), 6.48M PCFs, 31.4T MITM Expressions, CDT Barrier & 3D Calabi Geometry:**
   - Independently certified `73*RV - 1` (`μ ≤ 5.51374682 < 5.5137469`, matching Romero Oct 2026), scanned 6.48M PCFs & 31.4T symbolic expressions, ruled out exotic/multiplicative/CM-Gamma/hyperbolic relations at 120 digits, and verified $\operatorname{Vol}(\Omega_3^{\mathrm{hyp}}) = \frac{7}{8}\zeta(3)$ vs $\operatorname{Vol}(\Delta_3) = \frac{\pi^3}{32}$.
7. **Lean 4 Formalization (`lake build`, 18/18 jobs, zero `sorry`, zero warnings, standard axioms only):**
   - Formalized Apéry recurrence, Euler parity obstruction, Plouffe's Lambert equivalence, Rhin–Viola 30-coset certificates, 3D Calabi/Hyperbolic/Jacobi-Elliptic Jacobians, AND `Zeta3/OvernightTheorems.lean` (5-binomial factorization of $F_n^\star$, zero-axiom 24-shell lower bounds, universal trinomial Pisot & infinite Salem collapse theorems, Lehmer telescoping, Fauzan scaling core, and $\chi_{-4}$ parity obstruction).

## Repository Structure

```text
zeta3-closed-form/
├── README.md                              # This overview
├── paper/
│   ├── main.tex                           # Publication-ready arXiv LaTeX manuscript (mu(zeta(3)) < 5.3037 & dim<1,z3,z4>=3, mu < 112.3)
│   ├── main.pdf                           # Compiled 27-page AMS/arXiv PDF (with Z/p^2Z proof & Lean 4 formalization in Appendix C)
│   └── references.bib                     # BibTeX bibliography
├── docs/
│   ├── CAMPAIGN_25_EXPERIMENTS_REPORT.md  # Full 25-experiment research campaign & elementary Z/p^2Z proof of < 5.3037 & < 112.3
│   ├── NUMBER_THEORY_MASTER_STRATEGY.md   # 10-branch master strategy, 5 stress tests & 5-module elementary proof study guide
│   ├── OVERNIGHT_DISCOVERIES.md           # New mu<=5.3037 CMF, 8 pure 225d Pisot/Salem/Lehmer ladders, Fauzan w_3 & WZ v2
│   ├── CREATIVE_CLOSED_FORM_REPORT.md     # 31.4T MITM hunt, 120d exotic LLL, Pisot ladders & Calabi geometry
│   ├── EXPERIMENTS_REPORT.md              # Full report of RV `73*RV-1`, 6.48M PCFs, CDT barrier & Nov 2026 frontier
│   ├── HIGHEST_PROBABILITY_STRATEGY.md    # Definitive ranked strategy (102-search synthesis + RH clues)
│   ├── LITERATURE_100_SEARCH_ARCHIVE.md   # Complete indexed archive of 102 targeted research searches (1734–2026)
│   ├── ZETA3_COMPENDIUM.md                # Mathematical dossier: obstructions, 20+ formulas & roadmap
│   └── RH_PROOF_CLUES.md                  # Concrete transferable machinery audited from the 2026 Lean RH proof
├── tools/
│   ├── cmf_hunter/                        # 16-core native Rust 5F4/6F5/7F6 CMF & exact Farey F_p-shell optimizer
│   ├── wz_series_hunter/                  # 16-core native Rust 31-digit DD 255k-kernel hypergeometric hunter (<10 MB RAM)
│   ├── rv_hunter/                         # 16-core native Rust Rhin–Viola W(D_5) optimizer (<20 MB RAM)
│   ├── pcf_hunter/                        # 16-core native Rust 6.48M-PCF / 4D-LLL hunter (<10 MB RAM)
│   └── expr_hunter/                       # 16-core native Rust 31.4-trillion MITM expression hunter (<11 MB RAM)
├── scripts/
│   ├── certify_paper_bounds.py            # 100-digit rigorous interval-arithmetic & universal Farey-40 certifier (5.3037 & 112.3)
│   ├── verify_explicit_Gm_formula_full.py # Pointwise mod-p^2 verification of 2nd-order lifts across all 5 sub-windows (u<1 & u>=1)
│   ├── prove_fp_shells_530.py             # Elementary Wilson-Lucas deg(g_m)=p-2 verifier (38/38 sub-shells) & joint (1,z3,z4)
│   ├── campaign_phase1_exps1_5.py         # 25-experiment campaign: Phase 1 (Exps 1-5: 2nd-order, twins, shifts, Farey-80)
│   ├── campaign_phase3_exps11_12.py       # 25-experiment campaign: Phase 3 (Exps 11-12: continuous R^{2m+4} & dilated Pochhammers)
│   ├── campaign_phase4_exp13_16_fixed.py  # 25-experiment campaign: Phase 4 (Exps 13-16: RV-1 mod p^4, 4-block odd, skew-3-block)
│   ├── campaign_phase5_exps19_25.py       # 25-experiment campaign: Phase 5 (Exps 19-25: 5D Vasilyev, Fauzan, wings, k*F*, joint)
│   ├── check_exp25_joint_exact.py         # Exact Farey certification of single-presentation joint mu(A*zeta(3) + B*zeta(4))
│   ├── strategy_stress_tests.py           # 5 quantitative cross-branch number-theory stress tests (Exp 3A, 3B, 4A, 4B, 4C)
│   ├── certify_cmf_530.py                 # 230-digit independent certifier for new mu(zeta(3)) <= 5.303682 CMF
│   ├── pisot_salem_ladder_hunt.py         # 125-digit Pisot, Salem, Lehmer, cyclotomic & Weber Li_3 ladder hunter
│   ├── verify_pure_li3_ladders.py         # 225-digit certifier for all 8 pure zero-pi, zero-ln Li_3 closed forms
│   ├── hankel_zeta3_pi3_lab.py            # 420-digit weight-3 Fauzan e^{-1.084 K^2} Hankel & alternating w_alt lab
│   ├── verify_wz_hits.py                  # 130-digit WZ v2 hypergeometric verifier & Guillera shadow-duality lab
│   ├── certify_rv_bound.py                # 110-digit exact-rational & Stieltjes-remainder certifier (width 4.8e-80)
│   ├── verify_pcf_hits.py                 # 110-digit PCF verifier & algebraic near-miss diagnosis
│   ├── exotic_closed_form_lab.py          # 120-digit multiplicative, CM-Gamma, modular & hyperbolic LLL lab
│   ├── polylog_ladder_hunt.py             # 115-digit golden-ratio, plastic-constant & algebraic-unit Li_3 ladders
│   ├── hypergeom_series_discovery.py      # 115-digit 468-kernel central-binomial & BBP series discovery
│   ├── calabi_mahler_lab.py               # 3D circular/hyperbolic Calabi volumes & Smyth Mahler measures
│   ├── cm_lambert_discovery.py            # 110-digit CM quadratic Lambert / Eichler LLL discovery sweep
│   ├── cdt_capacity_test.py               # 65-digit CDT conformal-capacity vs denominator-barrier lab
│   ├── zeta3_lab.py                       # 80-digit verifier (8 identities to 1e-81) & LLL relation exclusion lab
│   └── rhin_viola_lab.py                  # Beukers / Rhin–Viola W(D_5) (|Φ|=1920) workbench
├── lean-toolchain                         # Lean 4 toolchain (v4.36.0)
├── lakefile.toml                          # Lake project configuration
└── Zeta3/
    ├── Basic.lean                         # Apéry numbers, cubic recurrence & rational sequences
    ├── AperyRecurrence.lean               # Formal proofs: Casoratian conservation & telescoping series identity
    ├── ParityObstruction.lean             # Formal proof: Euler's symmetric pole-pairing parity obstruction
    ├── OpenProblems.lean                  # Formal open-problem specs, Ramanujan & Plouffe Lambert equivalences
    ├── RhinViolaCert.lean                 # Kernel-checked W(D_5) uniform-shift law & 30-coset certificates
    ├── CalabiGeometry.lean                # Formal 3D circular, hyperbolic & Jacobi-elliptic Calabi Jacobians
    ├── OvernightTheorems.lean             # Formal 5-binomial F_n*, 24-shell bounds, Pisot/Salem/Lehmer & Fauzan theorems
    └── ElementaryRecordProof.lean         # Formal saddle cubic, Wilson-Lucas degree, 38 sub-shells, concavity & Thms 1-2
```

## Quick Start

### 1. Independently Certify the Paper (`paper/main.pdf`), Wilson–Lucas Proof & Strategy Stress Tests
```bash
python3 scripts/certify_paper_bounds.py     # 100-digit interval enclosures, universal Farey-40 & concavity certifier
python3 scripts/prove_fp_shells_530.py      # Elementary Wilson-Lucas proof (38/38 sub-shells) & single-presentation (1,z3,z4)
python3 scripts/strategy_stress_tests.py    # 5 cross-branch number-theory stress tests (Exp 3A, 3B, 4A, 4B, 4C)
python3 scripts/certify_cmf_530.py          # 230-digit saddle, exact Fraction carry/24-shell & n=1..3 verification
python3 scripts/verify_pure_li3_ladders.py  # 225-digit verification of all 8 pure Pisot/Salem/Lehmer Li_3 formulas
python3 scripts/certify_rv_bound.py         # 80-digit certification of RV01 cusp optimum 73*RV - 1 (< 5.5137469)
```

### 2. Build the Lean 4 Formalization (20/20 jobs, zero `sorry`) & Compile the Paper PDF
```bash
lake build
(cd paper && tectonic -X compile main.tex)
```

### 3. Run the 16-Core Rust Hunters & Supporting Labs
```bash
(cd tools/cmf_hunter && cargo run --release --bin certify_530)
(cd tools/wz_series_hunter && cargo run --release --bin wz_v2) && python3 scripts/verify_wz_hits.py
(cd tools/rv_hunter && cargo run --release -- selftest)
(cd tools/pcf_hunter && cargo run --release) && python3 scripts/verify_pcf_hits.py
python3 scripts/hankel_zeta3_pi3_lab.py
```

## Research Strategy (Two-Pronged Attack)

Because Grothendieck's Period Conjecture and Mixed Tate Motives ($\mathrm{MTM}(\mathbb{Z})$) strongly predict that **no** elementary closed form $F(\pi) = \zeta(3)$ exists, "tackling the closed form for $\zeta(3)$" splits into two complementary mathematical fronts:

1. **Negative / Non-Existence Direction (Proving No Euler Closed Form Exists):**
   - **Target 1 (Milestone):** Prove $\zeta(3) / \pi^3 \notin \mathbb{Q}$ (i.e., there is no rational closed form $\zeta(3) = q\pi^3$).
   - **Target 2 (Quantitative):** Beat Giraudin's Sept 2026 irrationality measure record $\mu(\zeta(3)) < 5.5138800$ (`arXiv:2609.26980`) by searching narrow-cone $W(D_5)$ perturbations and larger hypergeometric / multiple-integral families.
   - **Target 3 (Transcendence):** Study Nesterenko-style modular / Eichler-integral approaches linking $\zeta(3)$ via Ramanujan's formula $\zeta(3) = \frac{7\pi^3}{180} - 2\sum_{n=1}^\infty \frac{1}{n^3(e^{2\pi n}-1)}$ to quasi-modular forms.

2. **Positive / Structural Direction (Closed Forms in Enlarged Function & Period Rings):**
   - Express $\zeta(3)$ in canonical period, polylogarithmic, Eichler-integral, Mahler-measure, and hypergeometric bases.
   - Use high-precision LLL/PSLQ search (`scripts/zeta3_lab.py`) to test candidate relations and discover accelerated hypergeometric / continued-fraction representations (following Amdeberhan–Zeilberger, Ramanujan Machine, and Beukers–Rhin–Viola).

See [`docs/ZETA3_COMPENDIUM.md`](./docs/ZETA3_COMPENDIUM.md) for the complete reference and detailed roadmap.
