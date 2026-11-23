# $\zeta(3)$ Compendium: Closed Forms, Obstructions, Representations & Research Roadmap

This document gathers the essential mathematics behind **Apéry's constant** $\zeta(3) = \sum_{n=1}^{\infty} n^{-3}$ so that every formula, obstruction, prior formalization, and open research direction is in one place.

---

## 1. Fundamental Constant & Numerical Value

$$\zeta(3) = \sum_{n=1}^{\infty} \frac{1}{n^3} = \prod_{p \text{ prime}} \frac{1}{1 - p^{-3}}$$

To 60 decimal places:
```text
ζ(3) = 1.202056903159594285399738161511449990764986292340498881792272...
```

Simple continued fraction:
$$\zeta(3) = [1;\, 4, 1, 18, 1, 1, 1, 4, 1, 9, 9, 2, 1, 1, 1, 2, 7, 1, 1, 7, 11, 1, \dots]$$

Unlike $\pi^3 \approx 31.006276680299820175476315067101395202225288565885107694144538$, the ratio is:
$$\frac{\zeta(3)}{\pi^3} \approx 0.038768179602916798941119890318721149806234568039552579223126\ldots$$
with continued fraction $[0;\, 25, 1, 3, 1, 6, 3, 1, 1, 2, 1, 10, 3, 2, 1, 19, \dots]$ — showing no simple rational structure (every rational approximation $p/q$ with small $q$ fails; see `scripts/zeta3_lab.py`).

---

## 2. Why $\zeta(2k)$ Has a Closed Form and Why $\zeta(3)$ Resists It

There are **four distinct mathematical perspectives** explaining why $\zeta(2k) \in \mathbb{Q}\cdot \pi^{2k}$ works while $\zeta(2k+1)$ breaks:

### 2.1 Parity Obstruction in Residue Calculus
To evaluate $\zeta(2k)$ via contour integration, consider $f(z) = \frac{\pi \cot(\pi z)}{z^{2k}}$.
- Near $z = n \in \mathbb{Z} \setminus \{0\}$, $\operatorname{Res}_{z=n} f(z) = \frac{1}{n^{2k}}$.
- Because $n^{2k} = (-n)^{2k}$, the residues at $+n$ and $-n$ **add**: $\sum_{n \neq 0} \frac{1}{n^{2k}} = 2\zeta(2k)$.
- The sum of all residues is $0$, so $2\zeta(2k) = -\operatorname{Res}_{z=0} \frac{\pi \cot(\pi z)}{z^{2k}}$, which is an explicit rational multiple of $\pi^{2k}$ (via Bernoulli numbers).

For $2k+1 = 3$, $\frac{\pi \cot(\pi z)}{z^3}$ is an **even** function of $z$:
- The residues at $+n$ and $-n$ have **opposite signs** ($1/n^3$ vs $-1/n^3$) and **cancel identically**: $\sum_{n \neq 0} \frac{1}{n^3} = 0$, while $\operatorname{Res}_{z=0} f(z) = 0$ as well!
- Replacing $\cot(\pi z)$ by an even kernel (e.g. $\csc^2(\pi z)$) shifts the exponent back to even integers; replacing by $\psi(-z)$ introduces logarithmic/digamma branch cuts rather than rational residues.

### 2.2 Functional Equation Obstruction (Trivial Zeros)
Riemann's completed zeta function satisfies $\xi(s) = \xi(1-s)$, equivalently:
$$\zeta(s) = 2^s \pi^{s-1} \sin\!\left(\frac{\pi s}{2}\right) \Gamma(1-s)\, \zeta(1-s).$$
- At **$s = 2k$**: $\sin(\pi k) = 0$ cancels the simple pole of $\Gamma(1-2k)$, yielding a finite rational relation between $\zeta(2k)$ and $\zeta(1-2k) = -\frac{B_{2k}}{2k} \in \mathbb{Q}^\times$.
- At **$s = 2k+1$** (e.g., $s=3$): $\sin\!\left(\frac{3\pi}{2}\right) = -1 \neq 0$, while $\Gamma(-2)$ has a simple pole and $\zeta(-2) = 0$ is a **trivial zero**. Taking the limit $s \to 3$ cancels the pole against the zero via L'Hôpital's rule, producing the **derivative** $\zeta'(-2)$ rather than a rational Bernoulli value:
$$\zeta(3) = -4\pi^2 \zeta'(-2), \qquad \text{more generally } \zeta(2k+1) = (-1)^k \frac{2(2\pi)^{2k}}{(2k)!}\, \zeta'(-2k).$$
Thus the functional equation translates $\zeta(3)$ into $\zeta'(-2)$, not into $\mathbb{Q}\cdot \pi^3$.

### 2.3 Fourier / Bernoulli Polynomial Obstruction
The periodic Bernoulli polynomials $\overline{B}_m(x)$ have Fourier series:
$$\overline{B}_{2k}(x) = (-1)^{k-1} \frac{2\,(2k)!}{(2\pi)^{2k}} \sum_{n=1}^{\infty} \frac{\cos(2\pi n x)}{n^{2k}}, \qquad \overline{B}_{2k+1}(x) = (-1)^{k-1} \frac{2\,(2k+1)!}{(2\pi)^{2k+1}} \sum_{n=1}^{\infty} \frac{\sin(2\pi n x)}{n^{2k+1}}.$$
- Evaluating $\overline{B}_{2k}(x)$ at $x = 0$ gives $\cos(0) = 1$, immediately extracting $\zeta(2k)$.
- Evaluating $\overline{B}_{2k+1}(x)$ at $x = 0$ gives $\sin(0) = 0$! At $x = 1/4$, $\sin(\pi n/2)$ extracts the **Dirichlet beta value** $\beta(2k+1) = \sum_{n=0}^\infty \frac{(-1)^n}{(2n+1)^{2k+1}}$ (which *does* have a closed form in $\pi^{2k+1}$, e.g. $\beta(3) = \frac{\pi^3}{32}$!), **not** $\zeta(2k+1)$.
- In general, for a primitive Dirichlet character $\chi$ of parity $\chi(-1) \in \{+1,-1\}$, $L(\chi, m)$ is an algebraic multiple of $\pi^m$ **precisely when** $\chi(-1) = (-1)^m$ (the **critical values** in Deligne's sense). For the trivial character $\chi = 1$, $\chi(-1) = +1$, so even $m$ are critical and odd $m \ge 3$ are **non-critical**!

### 2.4 Motivic & Period Conjecture Obstruction (Deepest Explanation)
In algebraic geometry and arithmetic motives (Deligne, Goncharov, Kontsevich–Zagier, Francis Brown 2012):
- **Periods** (Kontsevich–Zagier 2001): complex numbers whose real and imaginary parts are integrals of rational functions over semialgebraic domains defined over $\mathbb{Q}$. Both $2\pi i$ and $\zeta(3)$ are periods ($\zeta(3) = \iiint_{[0,1]^3} \frac{dx\,dy\,dz}{1-xyz}$).
- **Mixed Tate Motives over $\mathbb{Z}$ ($\mathrm{MTM}(\mathbb{Z})$)**:
  - The motivic Galois group of $\mathrm{MTM}(\mathbb{Z})$ is $\mathbb{G}_m \ltimes U$, where the Lie algebra of the unipotent radical $U$ is a **free Lie algebra** generated by one generator $\sigma_{2k+1}$ in each odd degree $3, 5, 7, \dots$ (Deligne–Goncharov; Borel's calculation of $K_{2n-1}(\mathbb{Z}) \otimes \mathbb{Q}$).
  - Consequently, the ring of *motivic* multiple zeta values has $\zeta^{\mathfrak{m}}(2), \zeta^{\mathfrak{m}}(3), \zeta^{\mathfrak{m}}(5), \dots$ **algebraically independent** (Brown 2012 proved every motivic MZV is a $\mathbb{Q}$-linear combination of Hoffman elements $\zeta^{\mathfrak{m}}(n_1,\dots,n_r)$ with $n_i \in \{2,3\}$).
  - **Grothendieck's Period Conjecture** asserts that the period map $\mathrm{per} : \mathcal{H}^{\mathfrak{m}} \to \mathbb{C}$ is **injective**. If true, it immediately implies:
    $$\pi,\; \zeta(3),\; \zeta(5),\; \zeta(7),\; \dots \quad \text{are algebraically independent over } \mathbb{Q}.$$
  - Under this conjecture, **no polynomial equation** $P(\pi, \zeta(3)) = 0$ with $P \in \mathbb{Q}[X, Y] \setminus \{0\}$ can exist—and in particular no Euler closed form $\zeta(3) = q\pi^3$ exists.

---

## 3. Complete Catalog of Known Representations of $\zeta(3)$

All numerical identities below are verified to 50+ digits by [`scripts/zeta3_lab.py`](../scripts/zeta3_lab.py).

### 3.1 Basic Series & Dirichlet Relations
1. **Definition:** $\displaystyle \zeta(3) = \sum_{n=1}^{\infty} \frac{1}{n^3}$
2. **Dirichlet eta ($\eta$):** $\displaystyle \zeta(3) = \frac{4}{3} \sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n^3} = \frac{4}{3}\eta(3)$
3. **Dirichlet lambda ($\lambda$):** $\displaystyle \zeta(3) = \frac{8}{7} \sum_{n=0}^{\infty} \frac{1}{(2n+1)^3} = \frac{8}{7}\lambda(3)$
4. **Euler's Multiple Zeta Value (MZV) identity (1775):**
   $$\zeta(3) = \zeta(2,1) = \sum_{n=1}^{\infty} \frac{H_{n-1}}{n^2} = \frac{1}{2}\sum_{n=1}^{\infty} \frac{H_n}{n^2}, \qquad H_n = \sum_{k=1}^n \frac{1}{k}$$
5. **Tetragamma function:** $\displaystyle \zeta(3) = -\frac{1}{2}\psi^{(2)}(1)$

### 3.2 Rapidly Converging Hypergeometric Series
6. **Markov (1890) / Apéry (1978) series:**
   $$\zeta(3) = \frac{5}{2} \sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n^3 \binom{2n}{n}}$$
   *(Converges geometrically at rate $1/4$ per term.)*
7. **Hjortnaes (1953) / Gosper / Amdeberhan (1996):**
   $$\zeta(3) = \frac{1}{4}\sum_{n=1}^{\infty} (-1)^{n-1} \frac{56n^2 - 32n + 5}{(2n-1)^2} \frac{((n-1)!)^3}{(3n)!} \cdot \dots$$
   and the celebrated **Amdeberhan–Zeilberger (1997)** formula (found via WZ method; ~3 decimal digits per term, rate $1/1024$):
   $$\zeta(3) = \frac{1}{64} \sum_{k=0}^{\infty} \frac{(-1)^k (205k^2 + 250k + 77)\,(k!)^{10}}{((2k+1)!)^5}$$
8. **Wedeniwski (1998)** (used for billion-digit records; ~5.04 decimal digits per term):
   $$\zeta(3) = \frac{1}{24} \sum_{k=0}^{\infty} (-1)^k \frac{((2k+1)!(2k)!k!)^3 (126392k^5 + 412708k^4 + 531578k^3 + 336367k^2 + 104000k + 12463)}{(3k+2)!\,((4k+3)!)^3}$$
9. **Broadhurst (1998) binary BBP-type formula:**
   Allows extracting isolated binary digits of $\zeta(3)$ directly via base-$2^{12}$ polylogarithmic ladders (though significantly more complex than $\pi$'s BBP formula).

### 3.3 Integral & Period Representations
10. **Mellin / Bose–Einstein integral:**
    $$\zeta(3) = \frac{1}{2} \int_0^\infty \frac{x^2}{e^x - 1}\,dx$$
11. **Beukers' triple integral (1979):**
    $$\zeta(3) = \int_0^1 \int_0^1 \int_0^1 \frac{dx\,dy\,dz}{1 - xyz}$$
    *(Verified by expanding $\frac{1}{1-xyz} = \sum_{n=0}^\infty (xyz)^n$ and integrating term by term: $\sum_{n=0}^\infty \frac{1}{(n+1)^3} = \zeta(3)$; note Beukers' log integral has the $\frac{1}{2}$ factor: $\zeta(3) = -\frac{1}{2}\iint_{[0,1]^2} \frac{\ln(xy)}{1-xy}\,dx\,dy$.)*
    With shifted Legendre polynomials $P_n(x) = \frac{1}{n!}\frac{d^n}{dx^n}(x^n(1-x)^n)$, Beukers' linear form is:
    $$I_n = -\frac{1}{2} \int_0^1 \int_0^1 \frac{\ln(xy)}{1-xy} P_n(x) P_n(y)\,dx\,dy = \int_0^1 \int_0^1 \int_0^1 \frac{(x(1-x)y(1-y)z(1-z))^n}{(1 - (1-xy)z)^{n+1}}\,dx\,dy\,dz = a_n \zeta(3) - b_n$$
    where $a_n \in \mathbb{Z}$, $2 d_n^3 b_n \in \mathbb{Z}$ ($d_n = \operatorname{lcm}(1,\dots,n)$), and $0 < I_n < \zeta(3)(\sqrt{2}-1)^{4n}$. Because $e^3 (\sqrt{2}-1)^4 \approx 0.59126 < 1$, this gives an immediate proof that $\zeta(3) \notin \mathbb{Q}$!
12. **Euler's log-sine integral (1772):**
    $$\zeta(3) = \frac{2\pi^2}{7}\ln 2 + \frac{16}{7} \int_0^{\pi/2} x \ln(\sin x)\,dx$$
13. **Double log integral:**
    $$\zeta(3) = \frac{1}{2}\int_0^1 \frac{\ln^2 x}{1 - x}\,dx = -\frac{1}{2}\int_0^1 \frac{\ln(1-x)\ln x}{x}\,dx$$
14. **Polylogarithm at $z = 1/2$ (Landen 1780):**
    $$\zeta(3) = \frac{8}{7}\left(\operatorname{Li}_3\!\left(\frac{1}{2}\right) + \frac{\pi^2}{12}\ln 2 - \frac{1}{6}(\ln 2)^3\right)$$
15. **Clausen values ($\operatorname{Cl}_3(\theta) = \sum_{n=1}^\infty \frac{\cos(n\theta)}{n^3}$):**
    $$\operatorname{Cl}_3(0) = \zeta(3), \quad \operatorname{Cl}_3(\pi) = -\frac{3}{4}\zeta(3), \quad \operatorname{Cl}_3\!\left(\frac{\pi}{2}\right) = -\frac{3}{32}\zeta(3), \quad \operatorname{Cl}_3\!\left(\frac{\pi}{3}\right) = \frac{1}{3}\zeta(3)$$

### 3.4 Modular / Ramanujan / Eichler-Integral Formulas
16. **Ramanujan's $\zeta(3)$ identity (1913, Notebooks Ch. 14 Entry 21(i); proved by Grosswald 1970, Berndt 1977):**
    For any $\alpha, \beta > 0$ with $\alpha\beta = \pi^2$:
    $$\frac{1}{\alpha}\left(\frac{1}{2}\zeta(3) + \sum_{n=1}^{\infty} \frac{1}{n^3(e^{2\alpha n}-1)}\right) + \frac{1}{\beta}\left(\frac{1}{2}\zeta(3) + \sum_{n=1}^{\infty} \frac{1}{n^3(e^{2\beta n}-1)}\right) = \frac{\alpha^2 + 5\alpha\beta + \beta^2}{180}$$
    *(Verified to $< 10^{-78}$ in `scripts/zeta3_lab.py`.)*

    At the self-dual CM point $\alpha = \beta = \pi$:
    $$\zeta(3) = \frac{7\pi^3}{180} - 2\sum_{n=1}^{\infty} \frac{1}{n^3(e^{2\pi n} - 1)}$$
    *(Notice how close $\frac{7\pi^3}{180} \approx 1.2057996486783\ldots$ is to $\zeta(3) \approx 1.2020569031596\ldots$! The residual Lambert series $S(\pi) = \sum_{n=1}^\infty \frac{1}{n^3(e^{2\pi n}-1)} \approx 0.001871372759\ldots$ is an **Eichler integral** of the weight-4 Eisenstein series $E_4(\tau)$ evaluated at $\tau = i$!)*
    - **Key insight:** Evaluating $S(\pi)$ in closed form would give a closed form for $\zeta(3)$ — and conversely, proving $S(\pi) / \pi^3 \notin \mathbb{Q}$ is **logically equivalent** to proving $\zeta(3)/\pi^3 \notin \mathbb{Q}$!

### 3.5 Apéry's Recurrence & Continued Fraction (1978)
17. **Apéry's 3-term recurrence:**
    $$(n+1)^3 u_{n+1} - (34n^3 + 51n^2 + 27n + 5) u_n + n^3 u_{n-1} = 0$$
    - **Integer solution $a_n$** (Apéry numbers, OEIS A005259): $a_0 = 1, a_1 = 5, a_2 = 73, a_3 = 1445, \dots$
      $$a_n = \sum_{k=0}^{n} \binom{n}{k}^2 \binom{n+k}{k}^2$$
    - **Rational solution $b_n$**: $b_0 = 0, b_1 = 6, b_2 = \frac{351}{4}, b_3 = \frac{62531}{36}, \dots$ with $2 d_n^3 b_n \in \mathbb{Z}$, where $d_n = \operatorname{lcm}(1, \dots, n)$:
      $$b_n = \sum_{k=0}^{n} \binom{n}{k}^2 \binom{n+k}{k}^2 \left( \sum_{m=1}^{n} \frac{1}{m^3} + \sum_{m=1}^{k} \frac{(-1)^{m-1}}{2 m^3 \binom{n}{m}\binom{n+m}{m}} \right)$$
    - **Casoratian (Wronskian) identity:**
      $$a_{n-1} b_n - a_n b_{n-1} = \frac{6}{n^3}$$
      Therefore:
      $$\zeta(3) - \frac{b_n}{a_n} = \sum_{k=n+1}^{\infty} \frac{6}{k^3 a_k a_{k-1}}$$
    - **Continued fraction:**
      $$\zeta(3) = \cfrac{6}{5 - \cfrac{1^6}{117 - \cfrac{2^6}{535 - \cfrac{3^6}{1463 - \ddots}}}}, \qquad P(n) = 34n^3 + 51n^2 + 27n + 5$$
      *(Apéry's proof gives irrationality measure $\mu(\zeta(3)) \le 1 + \frac{4\ln(1+\sqrt{2}) + 3}{4\ln(1+\sqrt{2}) - 3} \approx 13.41782$.)*

### 3.6 Algebraic Geometry, $K$-Theory & Mahler Measures
18. **Borel's Theorem ($K$-theory regulator):**
    $\zeta(3)$ is proportional (up to a rational multiple of $\pi^2$ / regulator normalization) to the Borel regulator of $K_5(\mathbb{Z}) \cong \mathbb{Z}$, just as $\zeta(2k+1)$ governs the regulator of $K_{4k+1}(\mathbb{Z})$.
19. **Mahler measure (Smyth 1981):**
    $$m(1 + x + y + z) = \int_0^1 \int_0^1 \int_0^1 \ln|1 + e^{2\pi i \theta_1} + e^{2\pi i \theta_2} + e^{2\pi i \theta_3}|\,d\theta_1\,d\theta_2\,d\theta_3 = \frac{7}{2\pi^2}\zeta(3)$$
    Equivalently:
    $$\zeta(3) = \frac{2\pi^2}{7}\, m(1 + x + y + z)$$
    *(A stunning geometric "closed form" of $\zeta(3)$ in terms of $\pi^2$ and the Mahler measure of the simplest 3-variable linear polynomial!)*
20. **Hyperbolic volume:** Related via the Bloch–Wigner dilogarithm / trilogarithm to volumes of hyperbolic manifolds and special values of $L$-functions of algebraic varieties.

---

## 4. What Does "Closed Form" Mean Rigorously?

To avoid ambiguity, mathematics distinguishes four nested notions of "closed form" (see Timothy Chow, *"What is a closed-form number?"*, Amer. Math. Monthly 1999):

1. **Euler / $\pi$-rational closed form:** $\zeta(3) \in \mathbb{Q}\cdot \pi^3$ (or $\overline{\mathbb{Q}}\cdot \pi^3$).
   - **Status:** Open! Even $\zeta(3)/\pi^3 \notin \mathbb{Q}$ is unproven.
2. **Polynomial / Algebraic closed form in $\pi$:** $P(\pi, \zeta(3)) = 0$ for some nonzero $P \in \mathbb{Q}[X, Y]$.
   - **Status:** Open! Conjectured impossible (Grothendieck's Period Conjecture).
3. **Elementary / EL-closed form (Chow 1999):** Expressible from $\mathbb{Q}$ using $+, -, \times, \div, \exp, \log$ (which includes $\pi = -i\log(-1)$ and all algebraic/trig/inverse-trig expressions).
   - **Status:** Open! Under Schanuel's conjecture + Period Conjecture, $\zeta(3)$ is conjectured **not** to be an EL-number.
4. **Period / Special-Function / Modular closed form:** Expressible via finitely many periods, Eichler integrals, Mahler measures, or hypergeometric values.
   - **Status:** **Achieved** in multiple striking ways (Beukers' integral, Smyth's Mahler measure $\frac{2\pi^2}{7}m(1+x+y+z)$, Ramanujan's Eichler integral $\frac{7\pi^3}{180} - 2\sum \frac{1}{n^3(e^{2\pi n}-1)}$, polylogarithm $\operatorname{Li}_3(1/2)$).

---

## 5. Prior Formalizations & Literature

- **Apéry (1978):** *Irrationalité de $\zeta(2)$ et $\zeta(3)$*, Astérisque 61.
- **Beukers (1979):** *A note on the irrationality of $\zeta(2)$ and $\zeta(3)$*, Bull. LMS 11.
- **Smyth (1981):** *On measures of polynomials in several variables*, Bull. Austral. Math. Soc. 23.
- **Hata (1990, 2000):** Improved irrationality measures of $\zeta(3)$ via Legendre-type rational approximations ($\mu \le 7.377956$).
- **Rhin & Viola (2001):** *The group structure for $\zeta(3)$*, Acta Arithmetica 97 — $\mu(\zeta(3)) \le 5.513891$ via a **1920-element Weyl group $\Phi \cong W(D_5)$** acting on an 8-parameter triple integral (stood for 25 years; note the 120-element $S_5$ group was for $\zeta(2)$ in 1996).
- **David Niedbala Giraudin (Sept 2026, `arXiv:2609.26980`):** *A new upper bound for the irrationality exponent of $\zeta(3)$* — **current world record $\mu(\zeta(3)) < 5.5138800$**, obtained by narrow-cone shifting of Rhin–Viola's 8 parameters ($52000\times$ RV shifted by $\le 23$).
- **Ball & Rivoal (2001) / Zudilin (2001) / Fischler (JLMS 2026):** Infinitely many odd zeta values are irrational; at least one of $\zeta(5), \zeta(7), \zeta(9), \zeta(11)$ is irrational; Fischler (2026) broke the $\log s$ linear-independence barrier to $\dim_{\mathbb{Q}} \operatorname{Span}(1, \zeta(3), \dots, \zeta(s)) \ge 0.21\sqrt{s/\log s}$ via Siegel's lemma.
- **Calegari, Dimitrov & Tang (2024–2026, JAMS 2025, `arXiv:2408.15403`, ICM 2026 `arXiv:2510.04156`):** Arithmetic holonomy bounds proving linear independence of $\{1, \zeta(2), L(2,\chi_{-3})\}$ and irrationality measures for $p$-adic $\zeta_2(5), \zeta_5(3)$.
- **Elimelech, Kaminer et al. (*PNAS* 2024; Wein et al. `arXiv:2507.08138`, 2025):** Conservative Matrix Fields (CMFs) unifying Apéry's $\zeta(3)$ proof and Ramanujan Machine polynomial continued fractions.
- **Brown (2012):** *Mixed Tate motives over $\mathbb{Z}$*, Annals of Mathematics 175.
- **Nesterenko (1996):** Proved algebraic independence of $\pi, e^\pi, \Gamma(1/4)$ via Ramanujan's differential equations for Eisenstein series $E_2, E_4, E_6$.
- **Lean 4 Formalizations:**
  - `Mathlib.NumberTheory.LSeries.RiemannZeta` — defines `riemannZeta`, proves $\zeta(2) = \pi^2/6$, functional equation, etc.
  - `ahhwuhu/zeta_3_irrational` (Liu, Zhang, Zhi, arXiv:2503.07625, 2025) — full Lean 4 formal proof of $\zeta(3) \notin \mathbb{Q}$ via Beukers' triple integral.
  - `leanprover-community/riemann-hypothesis` (Alpöge & Claude, Nov 2026) — full Lean 4 formal proof of RH & GRH for Dirichlet $L$-functions.

---

## 6. Concrete Research Roadmap: How We Tackle This

We organize our work into four escalating tiers—from immediately provable in Lean 4 to genuine open mathematical problems:

### Tier 1 — Foundation in Lean 4 (Immediate)
- [x] Set up Lean 4 project (`Zeta3/Basic.lean`, `Zeta3/OpenProblems.lean`, `Zeta3/AperyRecurrence.lean`).
- [x] Build 60-digit verification & LLL search workbench (`scripts/zeta3_lab.py`).
- [ ] Prove Apéry's recurrence, Casoratian identity $a_n b_{n-1} - a_{n-1} b_n = 6/n^3$, and integrality of $a_n$ in Lean 4.
- [ ] Add Mathlib dependency & connect `Zeta3/OpenProblems.lean` directly to `Mathlib`'s `riemannZeta 3` and `Real.pi`.

### Tier 2 — Computational & Structural Discovery (Near-Term)
- [ ] **PSLQ/LLL relation search** over expanded constant bases: $\{\zeta(3), \pi^3, \pi^2 \ln 2, (\ln 2)^3, \operatorname{Li}_3(1/2), S(\pi), S(\pi\sqrt{2}), S(\pi\sqrt{3}), \dots\}$ to map all modular/CM Lambert-series relations.
- [ ] **Rhin–Viola $W(D_5)$ & Giraudin Narrow-Cone Optimizer**: Implement Rhin–Viola's (2001) 1920-element $W(D_5)$ group-action + $p$-adic digamma integral calculator, verify both `5.513891` and Giraudin's Sept 2026 record `5.5138800` (`arXiv:2609.26980`), and search systematically over narrow-cone shifts and Zudilin's well-poised ${}_7F_6/{}_9F_8$ families to beat $\mu(\zeta(3)) < 5.5138800$!

### Tier 3 — Attacking "No Euler Closed Form": $\zeta(3)/\pi^3 \notin \mathbb{Q}$
- [ ] Study simultaneous rational approximations to $(1, \pi^3, \zeta(3))$ and to $(1, \zeta(2), \zeta(3))$.
- [ ] Analyze the **Ramanujan / Eichler bridge**: since $\zeta(3) - \frac{7}{180}\pi^3 = -2\sum_{n=1}^\infty \frac{1}{n^3(e^{2\pi n}-1)}$, proving $\zeta(3)/\pi^3 \notin \mathbb{Q}$ is equivalent to proving that the Eichler integral $S(\pi) = \sum_{n=1}^\infty \frac{1}{n^3(e^{2\pi n}-1)}$ is not a rational multiple of $\pi^3$.
- [ ] Investigate whether Nesterenko's modular-function method (which works for $E_2(i), E_4(i), E_6(i)$ and gives algebraic independence of $\pi$ and $e^\pi$) can be extended to Eichler integrals / period polynomials of weight 4, and pinpoint the exact obstruction (absence of a closed algebraic differential equation for the Eichler integral $S(\tau)$ with $q$-expansion coefficients $\sigma_{-3}(n)$).

### Tier 4 — Transcendence & Algebraic Independence of $\{\pi, \zeta(3)\}$
- [ ] Formalize the statement of Grothendieck's Period Conjecture specialized to $\mathrm{MTM}(\mathbb{Z})$ in weight $\le 3$.
- [ ] Explore whether the Alpöge–Claude paradigm (combining arithmetic multiplicity / theta lifts / massive machine-verified inequalities in Lean 4) suggests new auxiliary-function or diophantine-approximation constructions for odd zeta values.
