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
| Best irrationality measure $\mu(\zeta(3))$ | Classical refereed record: `< 5.513891` (Rhin–Viola 2001). **2026 Preprint Frontier:** `< 5.5138800` (Giraudin, Sept 2026, `arXiv:2609.26980`); `< 5.513747` (J. V. Romero, Oct 12, 2026 — independently rediscovered & 80-digit/Lean-4 certified in this repo as `73*RV - 1` `< 5.5137469`); **`< 5.4197`** via Conservative Matrix Fields (Leibtag–Barkan–Kalisch–Weinbaum–Shalyt–Kaminer, Oct 20, 2026, `arXiv:2610.27799`, which also claims $\mathbb{Q}$-linear independence of $1, \zeta(3), \zeta(4)$). Conjectured $\mu = 2$. |
| Does $\zeta(3)$ have a closed form in a richer class (modular/Eichler/polylog/period/CMF)? | **Partially Yes** — exact identities exist via Ramanujan's & Plouffe's Lambert/Eichler series, Clausen/polylogarithms, Beukers' triple integral, Conservative Matrix Fields, and Smyth's Mahler measure |

## Key Results (See [`docs/EXPERIMENTS_REPORT.md`](./docs/EXPERIMENTS_REPORT.md))

1. **Rhin–Viola (2001) Theorem 5.1 Cusp Optimum — $\mu(\zeta(3)) < 5.5137469$ (Independent Rediscovery & 80-Digit/Lean-4 Certification):**
   - Using our 16-core native Rust engine (`tools/rv_hunter`), we decomposed Giraudin's Sept 2026 record (`arXiv:2609.26980`) as a positive uniform shift $51997\,P_{\mathrm{RV}} + 50\cdot\mathbf{1}$ inside the unique 2D $W(D_5)$-commuting 5-tie plane, and independently located the exact opposite-sign cusp minimum (matching J. V. Romero's Oct 12, 2026 sharpening `≤ 5.5137468... < 5.513747`):
     $$P_{73} = 73\,P_{\mathrm{RV}} - \mathbf{1} = (1167, 1240, 1386, 1094, 875, 802, 656, 948)$$
   - Certified to **80 significant digits** (width $4.77 \times 10^{-80}$, `scripts/certify_rv_bound.py`), cross-verified in Giraudin's unmodified `verify_zeta3.c`, and combinatorially kernel-proved in `Zeta3/RhinViolaCert.lean`:
     $$\mu(\zeta(3)) \le 5.51374681766592246858789304016140360188838686966604\ldots < 5.5137469$$
     *(Note: Beyond Rhin–Viola's family, Leibtag et al. `arXiv:2610.27799` (Oct 20, 2026) push the CMF preprint frontier to `< 5.4197`; meanwhile A. Fauzan (Zenodo `22826419`) / A. R. Hyman (Isabelle AFP `Zeta_5_Irrational`) and Fischler–Lai–Sprang–Zudilin (`arXiv:2611.16156`) advance $\zeta(5)$ and odd-zeta counts via Hankel determinants.)*
2. **6.48-Million Polynomial Continued Fraction (PCF) Hunt (`tools/pcf_hunter` + `scripts/verify_pcf_hits.py`):**
   - Scanned 6,479,306 PCFs across 6 families in 1.8s against $\zeta(3)$, $\zeta(3)/\pi^3$, $\pi^3$, $\zeta(3)/\pi^2$, recovering all 3 exact primitive cubic PCFs for $\zeta(3)$ to $10^{-118}$ and confirming zero exact hits for $\zeta(3)/\pi^3$.
3. **CM Lambert Sweep & Quantitative CDT Conformal-Capacity Barrier (`scripts/cm_lambert_discovery.py`, `scripts/cdt_capacity_test.py`):**
   - Verified $t(i/\sqrt{6}) = \lambda(2i) = (\sqrt{2}-1)^4$ to 65+ digits, ruled out $\zeta(3)/\pi^3 = p/q$ for $|q| < 4.58 \times 10^{25}$, and quantified the exact $\Gamma(2)$ cusp deficit (`log(16) - 3 = -0.2274`) vs $\Gamma_0(6)^+$ surplus (`+0.5255`).
4. **Creative Closed-Form Hunt (`docs/CREATIVE_CLOSED_FORM_REPORT.md`):**
   - **31.4-Trillion Symbolic Expression MITM Sweep (`tools/expr_hunter`):** 16-core Rust search across `3.14 × 10¹³` depth-4 expressions over 57 constants (`0.10s`, `10.35 MB` RAM) — all 3 self-tests pass (`0.00e0`), zero exact elementary hits for $\zeta(3)$.
   - **120-Digit Exotic & Multiplicative LLL (`scripts/exotic_closed_form_lab.py`):** Rules out power-products (rational powers/roots), CM Gamma periods ($\Gamma(1/3), \Gamma(1/4)$), Dedekind $\eta$, hyperbolic volumes ($G, V_{\mathrm{tet}}$), Euler $\gamma$, Glaisher $A$, and $e^{\zeta(3)}$ up to $\|c\|_2 \ge 1.5\times 10^7$–$2.2\times 10^9$.
   - **Pure Golden-Ratio & Plastic-Constant Closed Forms (`scripts/polylog_ladder_hunt.py`, verified to `< 6e-122`):** Exact 6-term positive-power $\operatorname{Li}_3$ formulas for $\zeta(3)$ with **zero $\pi$ and zero logarithms** in the two smallest Pisot units — the golden ratio $\phi$ and the cube-root plastic constant $\rho$ ($\rho^3 = \rho + 1$).
   - **3D Calabi vs. Hyperbolic Calabi & Mahler Measures (`scripts/calabi_mahler_lab.py` + `Zeta3/CalabiGeometry.lean`):** Numerically verified $\operatorname{Vol}(\Omega_3^{\mathrm{hyp}}) = \frac{7}{8}\zeta(3)$ (`1.55e-15`), $\operatorname{Vol}(\Delta_3) = \frac{\pi^3}{32}$ (`5.55e-16`), and Smyth's $m(1+x+y+z) = \frac{7}{2\pi^2}\zeta(3)$ (`7.22e-16`), with Lean 4 formal proofs of the circular, hyperbolic, and Jacobi-elliptic 3D Jacobians + global inverse.
5. **Lean 4 Formalization (`lake build`, 16/16 jobs, zero `sorry`, standard axioms only):**
   - Formalized Plouffe's Lambert Equivalence Theorem (`noEulerClosedForm_iff_plouffe_lambert`, `plouffe_implies_ramanujan`), complete kernel-verified 30-coset splitting certificates (`Zeta3/RhinViolaCert.lean`), and 3D Calabi / Hyperbolic / Jacobi-Elliptic Jacobian & inverse theorems (`Zeta3/CalabiGeometry.lean`).

## Repository Structure

```text
zeta3-closed-form/
├── README.md                              # This overview
├── docs/
│   ├── CREATIVE_CLOSED_FORM_REPORT.md     # 31.4T MITM hunt, 120d exotic LLL, Pisot ladders & Calabi geometry
│   ├── EXPERIMENTS_REPORT.md              # Full report of RV `73*RV-1`, 6.48M PCFs, CDT barrier & Nov 2026 frontier
│   ├── HIGHEST_PROBABILITY_STRATEGY.md    # Definitive ranked strategy (102-search synthesis + RH clues)
│   ├── LITERATURE_100_SEARCH_ARCHIVE.md   # Complete indexed archive of 102 targeted research searches (1734–2026)
│   ├── ZETA3_COMPENDIUM.md                # Mathematical dossier: obstructions, 20+ formulas & roadmap
│   └── RH_PROOF_CLUES.md                  # Concrete transferable machinery audited from the 2026 Lean RH proof
├── tools/
│   ├── rv_hunter/                         # 16-core native Rust Rhin–Viola W(D_5) optimizer (<20 MB RAM)
│   ├── pcf_hunter/                        # 16-core native Rust 6.48M-PCF / 4D-LLL hunter (<10 MB RAM)
│   └── expr_hunter/                       # 16-core native Rust 31.4-trillion MITM expression hunter (<11 MB RAM)
├── scripts/
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
    └── CalabiGeometry.lean                # Formal 3D circular, hyperbolic & Jacobi-elliptic Calabi Jacobians
```

## Quick Start

### 1. Certify $\mu(\zeta(3)) < 5.5137469$ to 80 Digits (0.5s, Python stdlib only)
```bash
python3 scripts/certify_rv_bound.py
```

### 2. Build the Lean 4 Formalization (zero `sorry`)
```bash
lake build
```

### 3. Run the 16-Core Rust Hunters & Supporting Labs
```bash
(cd tools/rv_hunter && cargo run --release -- selftest)
(cd tools/pcf_hunter && cargo run --release) && python3 scripts/verify_pcf_hits.py
python3 scripts/cdt_capacity_test.py
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
