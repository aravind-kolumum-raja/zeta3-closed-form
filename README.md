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
| Best irrationality measure $\mu(\zeta(3))$ | $\mu(\zeta(3)) \le 5.513891$ (Rhin–Viola 2001; conjectured $\mu(\zeta(3)) = 2$) |
| Does $\zeta(3)$ have a closed form in a richer class (modular/Eichler/polylog/period)? | **Partially Yes** — exact identities exist via Ramanujan's modular series, Clausen/polylogarithms, Beukers' triple integral, and Borel's $K_5(\mathbb{Z})$ regulator |

## Repository Structure

```text
zeta3-closed-form/
├── README.md                    # This overview
├── docs/
│   ├── ZETA3_COMPENDIUM.md      # Complete mathematical dossier, obstructions & attack roadmap
│   └── RH_PROOF_CLUES.md        # Concrete clues & transferable machinery from the 2026 Lean RH proof
├── scripts/
│   └── zeta3_lab.py             # 60-digit arbitrary-precision identity verifier, Apéry lab & LLL relation finder
├── lean-toolchain               # Lean 4 toolchain (v4.36.0)
├── lakefile.toml                # Lake project configuration
└── Zeta3/
    ├── Basic.lean               # Core definitions, Apéry sequences (a_n, b_n), and exact algebraic recurrences
    ├── OpenProblems.lean        # Formal statements of the closed-form, independence & irrationality-measure targets
    └── AperyRecurrence.lean     # Verified theorems on Apéry numbers, recurrences & Casoratian identities
```

## Quick Start

### 1. Run the 60-Digit Numerical & LLL Exploration Lab
No external packages required (uses Python 3 standard library `decimal` & `fractions`):
```bash
python3 scripts/zeta3_lab.py
```

### 2. Build the Lean 4 Formalization
```bash
lake build
```

## Research Strategy (Two-Pronged Attack)

Because Grothendieck's Period Conjecture and Mixed Tate Motives ($\mathrm{MTM}(\mathbb{Z})$) strongly predict that **no** elementary closed form $F(\pi) = \zeta(3)$ exists, "tackling the closed form for $\zeta(3)$" splits into two complementary mathematical fronts:

1. **Negative / Non-Existence Direction (Proving No Euler Closed Form Exists):**
   - **Target 1 (Milestone):** Prove $\zeta(3) / \pi^3 \notin \mathbb{Q}$ (i.e., there is no rational closed form $\zeta(3) = q\pi^3$).
   - **Target 2 (Quantitative):** Improve Rhin–Viola's irrationality measure $\mu(\zeta(3)) \le 5.513891$ by searching larger hypergeometric / multiple-integral / group-action families.
   - **Target 3 (Transcendence):** Study Nesterenko-style modular / Eichler-integral approaches linking $\zeta(3)$ via Ramanujan's formula $\zeta(3) = \frac{7\pi^3}{180} - 2\sum_{n=1}^\infty \frac{1}{n^3(e^{2\pi n}-1)}$ to quasi-modular forms.

2. **Positive / Structural Direction (Closed Forms in Enlarged Function & Period Rings):**
   - Express $\zeta(3)$ in canonical period, polylogarithmic, Eichler-integral, Mahler-measure, and hypergeometric bases.
   - Use high-precision LLL/PSLQ search (`scripts/zeta3_lab.py`) to test candidate relations and discover accelerated hypergeometric / continued-fraction representations (following Amdeberhan–Zeilberger, Ramanujan Machine, and Beukers–Rhin–Viola).

See [`docs/ZETA3_COMPENDIUM.md`](./docs/ZETA3_COMPENDIUM.md) for the complete reference and detailed roadmap.
