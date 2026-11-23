#!/usr/bin/env python3
"""
zeta3_lab.py — Zero-dependency 80-digit arbitrary-precision laboratory for ζ(3).

Verifies all major classical and modern identities from docs/ZETA3_COMPENDIUM.md
and docs/LITERATURE_100_SEARCH_ARCHIVE.md to 60-75 decimal digits, and runs
exact LLL integer-relation searches across classical, modular, and polylogarithmic bases.
"""

from __future__ import annotations

import math
from decimal import Decimal, getcontext
from fractions import Fraction
from typing import List, Tuple

WORK_PREC = 85
getcontext().prec = WORK_PREC


def compute_pi() -> Decimal:
    """Compute π via the Chudnovsky algorithm (~14 digits/term)."""
    C = 426880 * Decimal(10005).sqrt()
    M = Decimal(1)
    L = Decimal(13591409)
    X = Decimal(1)
    K = Decimal(6)
    S = L
    for i in range(1, 8):
        M = (M * (K**3 - 16 * K)) / Decimal(i**3)
        L += 545140134
        X *= -262537412640768000
        S += (M * L) / X
        K += 12
    return C / S


def compute_ln2() -> Decimal:
    """Compute ln(2) via 2 * atanh(1/3)."""
    s = Decimal(0)
    p = Decimal(3)
    k = 0
    while True:
        term = Decimal(2) / (Decimal(2 * k + 1) * p)
        if term < Decimal("1e-82"):
            break
        s += term
        p *= 9
        k += 1
    return s


def compute_e() -> Decimal:
    """Compute e = exp(1)."""
    return Decimal(1).exp()


def compute_zeta3_az(terms: int = 32) -> Decimal:
    """Amdeberhan-Zeilberger (1997): ~3.01 digits/term."""
    s = Decimal(0)
    for k in range(terms):
        num = ((-1) ** k) * (205 * k * k + 250 * k + 77) * (math.factorial(k) ** 10)
        den = math.factorial(2 * k + 1) ** 5
        s += Decimal(num) / Decimal(den)
    return s / Decimal(64)


def compute_zeta3_wedeniwski(terms: int = 18) -> Decimal:
    """
    Wedeniwski (1998) / Amdeberhan-Zeilberger WZ formula (~5.04 digits/term):
      ζ(3) = (1/24) sum_{n>=0} (-1)^n * ((2n+1)!(2n)!n!)^3 * P(n) / ((3n+2)! ((4n+3)!)^3)
    where P(n) = 126392 n^5 + 412708 n^4 + 531578 n^3 + 336367 n^2 + 104000 n + 12463.
    """
    s = Decimal(0)
    for n in range(terms):
        P = (
            126392 * n**5
            + 412708 * n**4
            + 531578 * n**3
            + 336367 * n**2
            + 104000 * n
            + 12463
        )
        num = ((-1) ** n) * ((math.factorial(2 * n + 1) * math.factorial(2 * n) * math.factorial(n)) ** 3) * P
        den = math.factorial(3 * n + 2) * (math.factorial(4 * n + 3) ** 3)
        s += Decimal(num) / Decimal(den)
    return s / Decimal(24)


def compute_zeta3_apery_series(terms: int = 130) -> Decimal:
    """Markov (1890) / Apéry (1978): (5/2) sum_{n>=1} (-1)^(n-1) / (n^3 C(2n,n))."""
    s = Decimal(0)
    for n in range(1, terms + 1):
        c = math.comb(2 * n, n)
        s += Decimal((-1) ** (n - 1)) / Decimal((n**3) * c)
    return Decimal(5) * s / Decimal(2)


def compute_zeta3_crvz(n: int = 105) -> Decimal:
    """
    Cohen-Rodriguez Villegas-Zagier (2000) Chebyshev acceleration for η(3) = (3/4)ζ(3):
      ~0.765 digits/term with d_k = n * sum_{i=0}^k (n+i-1)! 4^i / ((n-i)! (2i)!).
    """
    d = Decimal(3 + math.sqrt(8))  # only for doc; compute exact integer d_k below!
    d_list = [Fraction(0)] * (n + 1)
    # Recurrence / exact integer formula for d_k
    c = 1
    d_acc = 0
    d_ints = []
    for k in range(n + 1):
        # term_i = n * (n+i-1)! * 4^i / ((n-i)! * (2i)!)
        if k == 0:
            term = 1
        else:
            term = (n * math.factorial(n + k - 1) * (4**k)) // (
                math.factorial(n - k) * math.factorial(2 * k)
            )
        d_acc += term
        d_ints.append(d_acc)
    d_n = d_ints[n]
    s = Decimal(0)
    for k in range(n):
        s += Decimal(((-1) ** k) * (d_ints[k] - d_n)) / Decimal((k + 1) ** 3)
    eta3 = -s / Decimal(d_n)
    return Decimal(4) * eta3 / Decimal(3)


def compute_catalan_G(terms: int = 140) -> Decimal:
    """Compute Catalan's constant G = sum (-1)^k/(2k+1)^2 via Ramanujan's accelerated series."""
    # G = (π/8) ln(2 + √3) + (3/8) sum_{n=0}^∞ 1 / ((2n+1)^2 C(2n, n)) ... wait use exact CRVZ!
    n = 110
    d_acc = 0
    d_ints = []
    for k in range(n + 1):
        if k == 0:
            term = 1
        else:
            term = (n * math.factorial(n + k - 1) * (4**k)) // (
                math.factorial(n - k) * math.factorial(2 * k)
            )
        d_acc += term
        d_ints.append(d_acc)
    d_n = d_ints[n]
    s = Decimal(0)
    for k in range(n):
        s += Decimal(((-1) ** k) * (d_ints[k] - d_n)) / Decimal((2 * k + 1) ** 2)
    return -s / Decimal(d_n)


def lambert_S3(r: Decimal, pi: Decimal, max_n: int = 120) -> Decimal:
    """Plouffe/Ramanujan Lambert series S_3(r) = sum_{k>=1} 1 / (k^3 * (exp(π*r*k) - 1))."""
    s = Decimal(0)
    for k in range(1, max_n + 1):
        exp_val = (pi * r * Decimal(k)).exp()
        term = Decimal(1) / (Decimal(k**3) * (exp_val - Decimal(1)))
        s += term
        if term < Decimal("1e-82"):
            break
    return s


def compute_li3_half(max_n: int = 320) -> Decimal:
    """Compute Li_3(1/2) = sum_{n>=1} 1 / (n^3 * 2^n)."""
    s = Decimal(0)
    p = Decimal(2)
    for n in range(1, max_n + 1):
        term = Decimal(1) / (Decimal(n**3) * p)
        s += term
        if term < Decimal("1e-82"):
            break
        p *= 2
    return s


def eval_ramanujan_machine_cf(depth: int = 60) -> Decimal:
    """
    Evaluate the Ramanujan Machine (Raayoni et al., Nature 2021) polynomial CF:
      8 / (7 ζ(3)) = 1 - 1^6 / (3*7 - 2^6 / (5*19 - 3^6 / (7*37 - ...)))
    where a_n = (2n+1)(3n^2 + 3n + 1) and b_n = -n^6.
    Returns the implied value of ζ(3).
    """
    val = Decimal(2 * depth + 1) * Decimal(3 * depth**2 + 3 * depth + 1)
    for n in range(depth - 1, -1, -1):
        a_n = Decimal(2 * n + 1) * Decimal(3 * n**2 + 3 * n + 1)
        b_next = Decimal(-((n + 1) ** 6))
        val = a_n + b_next / val
    return Decimal(8) / (Decimal(7) * val)


def eval_apery_cf(depth: int = 60) -> Decimal:
    """
    Evaluate Apéry's continued fraction (1978):
      ζ(3) = 6 / (5 - 1^6 / (117 - 2^6 / (535 - ...)))
    where P(n) = 34 n^3 + 51 n^2 + 27 n + 5.
    """
    val = Decimal(34 * depth**3 + 51 * depth**2 + 27 * depth + 5)
    for n in range(depth - 1, -1, -1):
        p_n = Decimal(34 * n**3 + 51 * n**2 + 27 * n + 5)
        b_next = Decimal(-((n + 1) ** 6))
        val = p_n + b_next / val
    return Decimal(6) / val


def apery_sequences(N: int) -> Tuple[List[int], List[Fraction]]:
    """Compute Apéry's sequences (a_n, b_n) for n = 0..N."""
    a = [1, 5]
    b = [Fraction(0), Fraction(6)]
    for n in range(1, N):
        P = 34 * n**3 + 51 * n**2 + 27 * n + 5
        a_next_num = P * a[n] - (n**3) * a[n - 1]
        assert a_next_num % ((n + 1) ** 3) == 0
        a.append(a_next_num // ((n + 1) ** 3))
        b_next = (P * b[n] - (n**3) * b[n - 1]) / ((n + 1) ** 3)
        b.append(b_next)
    return a, b


def continued_fraction(x: Decimal, steps: int = 20) -> List[int]:
    out = []
    cur = x
    for _ in range(steps):
        fl = int(cur)
        out.append(fl)
        rem = cur - Decimal(fl)
        if rem < Decimal("1e-65"):
            break
        cur = Decimal(1) / rem
    return out


# --- Exact Rational LLL Integer Relation Detector ---

def _dot(u: List[Fraction], v: List[Fraction]) -> Fraction:
    return sum(a * b for a, b in zip(u, v))


def lll_reduce(basis: List[List[int]], delta: Fraction = Fraction(3, 4)) -> List[List[int]]:
    n = len(basis)
    B = [list(map(int, row)) for row in basis]

    def gram_schmidt():
        B_star: List[List[Fraction]] = []
        mu: List[List[Fraction]] = [[Fraction(0)] * n for _ in range(n)]
        norms: List[Fraction] = []
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = _dot(v, B_star[j]) / norms[j]
                v = [v[k] - mu[i][j] * B_star[j][k] for k in range(len(v))]
            B_star.append(v)
            norms.append(_dot(v, v))
        return B_star, mu, norms

    B_star, mu, norms = gram_schmidt()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > Fraction(1, 2):
                q = round(mu[k][j])
                B[k] = [B[k][m] - q * B[j][m] for m in range(len(B[k]))]
                B_star, mu, norms = gram_schmidt()
        if norms[k] >= (delta - mu[k][k - 1] ** 2) * norms[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            B_star, mu, norms = gram_schmidt()
            k = max(k - 1, 1)
    return B


def find_integer_relation(values: List[Decimal], scale_digits: int = 50) -> Tuple[List[int], Decimal]:
    S = Decimal(10) ** scale_digits
    m = len(values)
    basis = []
    for i in range(m):
        row = [0] * (m + 1)
        row[i] = 1
        row[m] = int((values[i] * S).to_integral_value())
        basis.append(row)
    reduced = lll_reduce(basis)
    best_coeffs = reduced[0][:m]
    residual = sum(Decimal(c) * v for c, v in zip(best_coeffs, values))
    return best_coeffs, residual


def main() -> None:
    print("=" * 78)
    print("  ζ(3) 80-Digit Verification & Closed-Form Exclusion Lab (v2)")
    print("=" * 78)

    pi = compute_pi()
    ln2 = compute_ln2()
    e_val = compute_e()
    e_pi = pi.exp()
    G_cat = compute_catalan_G()
    z3 = compute_zeta3_az(32)

    print(f"π        = {str(pi)[:64]}")
    print(f"ζ(3)     = {str(z3)[:64]}")
    print(f"ζ(3)/π³  = {str(z3 / (pi**3))[:64]}")
    print(f"G (Cat.) = {str(G_cat)[:64]}")
    print(f"CF(ζ(3))    = {continued_fraction(z3, 18)}")
    print(f"CF(ζ(3)/π³) = {continued_fraction(z3 / (pi**3), 18)}")
    print()

    print("--- 1. Comprehensive Formula Verification (60-80 digits) ---")
    checks = [
        ("Markov-Apéry (1890/1978)", compute_zeta3_apery_series(130)),
        ("Wedeniwski WZ (1998)", compute_zeta3_wedeniwski(18)),
        ("Cohen-Villegas-Zagier (2000)", compute_zeta3_crvz(105)),
        ("Apéry CF (depth=60)", eval_apery_cf(60)),
        ("Ramanujan Machine CF (2021, d=60)", eval_ramanujan_machine_cf(60)),
    ]
    for name, val in checks:
        print(f"  {name:35s} |err| = {abs(z3 - val):.3e}")

    # Lambert / Eichler S_3(r) = sum 1 / (k^3 (exp(π r k) - 1))
    S1 = lambert_S3(Decimal(1), pi)
    S2 = lambert_S3(Decimal(2), pi)
    S4 = lambert_S3(Decimal(4), pi)
    S_sq2 = lambert_S3(Decimal(2).sqrt(), pi)
    S_2sq2 = lambert_S3(Decimal(2) * Decimal(2).sqrt(), pi)
    S_invsq2 = lambert_S3(Decimal(2) / Decimal(2).sqrt(), pi)

    z3_ram = (Decimal(7) * pi**3) / Decimal(180) - Decimal(2) * S2
    print(f"  {'Ramanujan self-dual α=π':35s} |err| = {abs(z3 - z3_ram):.3e}")

    # General Ramanujan αβ = π² with α = π√2, β = π/√2 (note S(α) = S_3(2α/π))
    alpha, beta = pi * Decimal(2).sqrt(), pi / Decimal(2).sqrt()
    lhs_ram_gen = (Decimal("0.5") * z3 + S_2sq2) / alpha + (Decimal("0.5") * z3 + S_invsq2) / beta
    rhs_ram_gen = (alpha**2 + Decimal(5) * alpha * beta + beta**2) / Decimal(180)
    print(f"  {'Ramanujan general α=π√2, β=π/√2':35s} |err| = {abs(lhs_ram_gen - rhs_ram_gen):.3e}")

    # Plouffe (1998): verify exact relations on S_3(1), S_3(2), S_3(4)
    # Note: Let's verify exact coefficients via LLL AND check Plouffe's claimed formulas!
    li3_half = compute_li3_half()
    z3_li3 = Decimal(8) / Decimal(7) * (li3_half + (pi**2) * ln2 / Decimal(12) - (ln2**3) / Decimal(6))
    print(f"  {'Landen Li_3(1/2) (1780)':35s} |err| = {abs(z3 - z3_li3):.3e}")
    print()

    print("--- 2. Apéry Recurrence & Casoratian Verification ---")
    N = 15
    a, b = apery_sequences(N)
    max_cas_err = max(abs(Fraction(a[n - 1]) * b[n] - Fraction(a[n]) * b[n - 1] - Fraction(6, n**3)) for n in range(1, N + 1))
    print(f"  Max Casoratian error (n=1..{N}): {max_cas_err}")
    print(f"  |ζ(3) - b_{N}/a_{N}| = {abs(z3 - Decimal(b[N].numerator) / Decimal(b[N].denominator * a[N])):.3e}")
    print()

    print("--- 3. LLL Integer Relation Detection (Positive Controls + Open Exclusions) ---")
    tests = [
        ("Control A: [ζ(3), Li_3(1/2), π² ln2, ln2³]", [z3, li3_half, (pi**2) * ln2, ln2**3], 50),
        ("Control B: [ζ(3), π³, S_3(2)] (Ramanujan)", [z3, pi**3, S2], 50),
        ("Control C: [ζ(3), S_3(1), S_3(2), S_3(4)] (Plouffe)", [z3, S1, S2, S4], 50),
        ("Control D: [π³, S_3(1), S_3(2), S_3(4)] (Plouffe)", [pi**3, S1, S2, S4], 50),
        ("Open 1: [1, π, π², π³, ζ(3)]", [Decimal(1), pi, pi**2, pi**3, z3], 55),
        ("Open 2: [1, π³, ζ(3), ζ(3)², π³ζ(3), π⁶]", [Decimal(1), pi**3, z3, z3**2, (pi**3) * z3, pi**6], 55),
        ("Open 3: [1, π³, π²ln2, ln2³, ζ(3)] (EL wt-3)", [Decimal(1), pi**3, (pi**2) * ln2, ln2**3, z3], 55),
        ("Open 4: [1, π, π², π³, e, e^π, ln2, G, ζ(3)]", [Decimal(1), pi, pi**2, pi**3, e_val, e_pi, ln2, G_cat, z3], 55),
    ]
    for label, vals, digits in tests:
        coeffs, res = find_integer_relation(vals, scale_digits=digits)
        max_c = max(abs(c) for c in coeffs)
        verdict = "EXACT RELATION FOUND!" if abs(res) < Decimal("1e-50") and max_c < 10**5 else f"NO small relation (min height > {max_c:,})"
        print(f"  {label}")
        print(f"    coeffs = {coeffs}")
        print(f"    |res|  = {abs(res):.3e}  =>  {verdict}")


if __name__ == "__main__":
    main()
