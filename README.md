# Improved Irrationality Measure $\mu(\zeta(3)) < 5.3037$ and Single-Presentation Linear Independence of $1, \zeta(3), \zeta(4)$ via Elementary Congruences

**Author:** Gemini 4 Neon — steered by Aravind Kolumum Raja ([Recommendium Labs](mailto:aravind@recommendium.com))  
**Manuscript:** [`paper/main.pdf`](./paper/main.pdf) ([LaTeX source](./paper/main.tex))  
**Status:** Open research preprint (not yet peer-reviewed; independent specialist review of Sections 3–4 is warmly welcomed)

---

## Main Results

Using new balanced parameters on the three-block ${}_5F_4(1)$ rational-function template of Leibtag et al. (`arXiv:2610.27799`, Appendix B.4), namely
$$F_n^{\star}(t) = \frac{(t - 13n)_{13n}\,(t - 15n)_{13n}\,(t - 17n)_{13n}}{(t + 2n)_{21n + 1}\,(t + 3n)_{18n + 1}},$$
we prove:

1. **Theorem 1.2 (Irrationality Measure of $\zeta(3)$):**
   $$\mu(\zeta(3)) \le 5.303681511 < 5.3037,$$
   improving upon the classical bound $5.513891$ of Rhin and Viola (2001) and the recent bound $5.4197$ of Leibtag, Barkan, Kalisch, Weinbaum, Shalyt, and Kaminer (`arXiv:2610.27799`, October 2026) by more than $0.1160$. *(First-order modulo-$p$ shells alone already yield $\mu(\zeta(3)) < 5.3193$. Given active unrefereed preprints/repositories in late 2026—including Davidson's November 2026 claim of $4.974$—we frame this bound strictly against the cited published/arXiv literature and emphasize the direct, self-contained proof method.)*

2. **Theorem 1.3 (Single-Presentation Linear Independence of $1, \zeta(3), \zeta(4)$):**
   The numbers $1$, $\zeta(3)$, and $\zeta(4) = \pi^4/90$ are linearly independent over $\mathbb{Q}$, with quantitative bound
   $$\mu\bigl(A\zeta(3) + B\zeta(4)\bigr) \le 112.248844 < 112.3 \qquad \forall\,(A, B) \in \mathbb{Q}^2 \setminus \{(0, 0)\},$$
   obtained from $F_n^{\star}(t)$ alone without a second dual presentation. *(First-order modulo-$p$ shells alone already yield $< 119.7$. Note that `arXiv:2610.27799` achieves a tighter joint bound $< 66.3$ by combining two distinct presentations.)*

### Why the Proof Is Direct and Self-Contained

The proof avoids Conservative Matrix Fields, Meijer $G$-functions, and multi-dimensional Mellin–Barnes integrals:

- **One-dimensional steepest descent (Section 3):** Both remainders $R_{3,n}$ and $R_{4,n}$ are represented by a single vertical Sommerfeld–Watson integral whose log-magnitude $\psi(y)$ is proved algebraically to be strictly concave on $(0, \infty)$ ($\psi''(y) < 0$ via a $9$-term negative polynomial in $w = y^2$).
- **Elementary modulo-$p$ shells (Section 4.2):** All $24$ first-order shells ($38$ uniform-carry sub-shells; $88$ open Farey-$40$ cells) follow universally for all primes $p > \max(40, \sqrt{40n})$ with $p \nmid n$ from an explicit Wilson–Lucas polynomial $\mathcal{R}_m(x) \in \mathbb{F}_p[x]$ of exact degree $p - 2 < p - 1$.
- **Elementary modulo-$p^2$ shells (Section 4.3):** Both second-order shells $I_{12} = (\frac{7}{15}, \frac{11}{23})$ and $I_{23} = (\frac{16}{17}, \frac{20}{21})$ ($5$ sub-windows; $14$ open Farey-$40$ cells) follow for all primes $p > \max(77, \sqrt{40n})$ with $p \nmid n$ via smooth modulo-$p^2$ binomial and simple-pole wing lifts (Wolstenholme's theorem), a $6$-line Faulhaber master congruence $\frac{1}{p}\sum_{x=0}^{p-1} F(x) \equiv \sum_{x=0}^{p-1} F(x)H_x \pmod p$ ($\deg F \le p - 2$), and an explicit degree-balanced quotient rational function $\widetilde{G}_q(t)$ with a triple zero at $t = 0$ resolving the $P_{4,n}$ bucket split on $I_{23}$.

---

## Independent Verification & Quick Start

All numerical enclosures, finite-field/modulo-$p^2$ identities, and auxiliary Lean 4 algebraic checks can be verified independently with zero external Python packages and zero external Lean libraries:

```bash
# 1. 100-digit interval-arithmetic certifier (zero external dependencies):
#    Certifies all 38 1st-order sub-shells (88 Farey-40 cells), all 14 Farey-40 cells of I_12 U I_23,
#    saddle roots, exact interval concavity, alternating-Bernoulli digamma enclosures, and final bounds.
python3 scripts/certify_paper_bounds.py

# 2. Direct mod-p^k verification of all 38 first-order sub-shells (u < 1 and periodic u >= 1):
python3 scripts/prove_fp_shells_530.py

# 3. Direct mod-p^4 pointwise verification of the second-order modulo-p^2 lifts (all sub-windows):
python3 scripts/verify_explicit_Gm_formula_full.py

# 4. 230-digit small-n verification (n = 1, 2, 3; Appendix B, Table 4):
python3 scripts/certify_cmf_530.py

# 5. Auxiliary Lean 4 algebraic & rational-arithmetic checks (zero sorry, zero extra axioms):
#    Note: checks auxiliary identities (saddle cubic, Delta periodicity, degree identity,
#    concavity polynomial at c=43/4, crossing identities, sub-shell defects, u<1 sums, and
#    cross-multiplied integer inequalities); does NOT formalize zeta(3) or the analytic/congruence proofs.
lake build

# 6. Recompile the LaTeX manuscript (requires Tectonic or TeX Live):
(cd paper && tectonic -X compile main.tex)
```

---

## Repository Structure

```text
zeta3-closed-form/
├── README.md                                 # Overview and reproduction instructions
├── paper/
│   ├── main.tex                              # AMS LaTeX manuscript
│   ├── main.pdf                              # Compiled PDF
│   └── references.bib                        # BibTeX bibliography
├── Zeta3.lean                                # Root Lean 4 module
├── Zeta3/
│   └── ElementaryRecordProof.lean            # Zero-dependency pure-ASCII Lean 4 checks (v4.36.0)
├── lean-toolchain                            # Lean toolchain pin
├── lakefile.toml                             # Lake build config
├── lake-manifest.json
├── scripts/
│   ├── certify_paper_bounds.py               # 100-digit interval-arithmetic & Farey-40 certifier
│   ├── prove_fp_shells_530.py                # First-order Wilson-Lucas mod-p^k verifier (38 sub-shells)
│   ├── verify_explicit_Gm_formula_full.py    # Second-order modulo-p^2 pointwise verifier
│   ├── certify_cmf_530.py                    # 230-digit exact small-n (n=1,2,3) & saddle verifier
│   ├── campaign_phase1_exps1_5.py            # Experimental audit: Farey-80 scan, twins, shifts (Exps 1-5)
│   ├── campaign_phase3_exps11_12.py          # Experimental audit: continuous R^{2m+4} & dilated blocks (Exps 11-12)
│   ├── campaign_phase4_exp13_16_fixed.py     # Experimental audit: RV-1 mod p^4 & 4-block odd (Exps 13-16)
│   ├── campaign_phase5_exps19_25.py          # Experimental audit: 5D, wings, k*F* tie-breaks (Exps 19-25)
│   └── check_exp25_joint_exact.py            # Exact Farey single-presentation joint checks (Exp 25)
├── tools/
│   └── cmf_hunter/                           # 16-core native Rust search & Farey certifier for F_n*(t)
└── docs/
    └── CAMPAIGN_25_EXPERIMENTS_REPORT.md     # Complete 25-experiment computational report
```

## Transparency & Contact

The mathematical proofs, search code, and companion checks were authored by **Gemini 4 Neon**, initiated and steered as an independent research project by **Aravind Kolumum Raja** (Recommendium Labs); see Section 6 of [`paper/main.pdf`](./paper/main.pdf) for the full division of labor and preprint status. Questions, referee comments, and independent verification reports are welcome at `aravind@recommendium.com`.
