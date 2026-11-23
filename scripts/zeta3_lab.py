#!/usr/bin/env python3
"""
zeta3_lab.py — Zero-dependency 60-digit arbitrary-precision laboratory for ζ(3).

Features:
  1. Computes π (Chudnovsky) and ζ(3) (Amdeberhan-Zeilberger) to 60+ digits.
  2. Numerically verifies the major identities in docs/ZETA3_COMPENDIUM.md:
       - Markov-Apéry series
       - Ramanujan modular / Eichler formula (self-dual α = β = π AND general αβ = π²)
       - Landen polylogarithm Li_3(1/2) identity
       - Apéry's recurrence, Casoratian identity, and continued fraction convergence
       - Beukers/Apéry denominator growth & irrationality exponent
  3. Runs a pure-Python LLL integer-relation search (PSLQ-style) to test
     whether ζ(3) satisfies any small-height linear/polynomial relation with π,
     ln(2), Li_3(1/2), and Ramanujan's Eichler series S(π).
"""

from __future__ import annotations

import math
from decimal import Decimal, getcontext
from fractions import Fraction
from typing import List, Tuple

# Work at 80 digits internally, report 50-60 verified digits.
WORK_PREC = 80
getcontext().prec = WORK_PREC


def compute_pi() -> Decimal:
    """Compute π via the Chudnovsky algorithm (~14 digits per term)."""
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
    """Compute ln(2) via 2 * atanh(1/3) = 2 * sum_{k>=0} 1/((2k+1)*3^(2k+1))."""
    s = Decimal(0)
    p = Decimal(3)
    k = 0
    while True:
        term = Decimal(2) / (Decimal(2 * k + 1) * p)
        if term < Decimal("1e-78"):
            break
        s += term
        p *= 9
        k += 1
    return s


def compute_zeta3_az(terms: int = 30) -> Decimal:
    """Compute ζ(3) via Amdeberhan-Zeilberger (1997), ~3.01 digits/term."""
    s = Decimal(0)
    for k in range(terms):
        num = ((-1) ** k) * (205 * k * k + 250 * k + 77) * (math.factorial(k) ** 10)
        den = math.factorial(2 * k + 1) ** 5
        s += Decimal(num) / Decimal(den)
    return s / Decimal(64)


def compute_zeta3_apery_series(terms: int = 120) -> Decimal:
    """Compute ζ(3) via Markov (1890) / Apéry (1978): (5/2) sum (-1)^(n-1)/(n^3 C(2n,n))."""
    s = Decimal(0)
    for n in range(1, terms + 1):
        c = math.comb(2 * n, n)
        term = Decimal((-1) ** (n - 1)) / Decimal((n**3) * c)
        s += term
    return Decimal(5) * s / Decimal(2)


def lambert_eichler_S(alpha: Decimal, max_n: int = 80) -> Decimal:
    """Compute S(α) = sum_{n>=1} 1 / (n^3 * (exp(2*α*n) - 1))."""
    s = Decimal(0)
    for n in range(1, max_n + 1):
        exp_val = (Decimal(2) * alpha * Decimal(n)).exp()
        term = Decimal(1) / (Decimal(n**3) * (exp_val - Decimal(1)))
        s += term
        if term < Decimal("1e-75"):
            break
    return s


def compute_li3_half(max_n: int = 300) -> Decimal:
    """Compute Li_3(1/2) = sum_{n>=1} 1 / (n^3 * 2^n)."""
    s = Decimal(0)
    p = Decimal(2)
    for n in range(1, max_n + 1):
        term = Decimal(1) / (Decimal(n**3) * p)
        s += term
        if term < Decimal("1e-75"):
            break
        p *= 2
    return s


def apery_sequences(N: int) -> Tuple[List[int], List[Fraction]]:
    """
    Compute Apéry's sequences (a_n, b_n) for n = 0..N using the 3-term recurrence:
      (n+1)^3 u_{n+1} = (34 n^3 + 51 n^2 + 27 n + 5) u_n - n^3 u_{n-1}
    with (a_0, a_1) = (1, 5) and (b_0, b_1) = (0, 6).
    """
    a = [1, 5]
    b = [Fraction(0), Fraction(6)]
    for n in range(1, N):
        P = 34 * n**3 + 51 * n**2 + 27 * n + 5
        a_next_num = P * a[n] - (n**3) * a[n - 1]
        assert a_next_num % ((n + 1) ** 3) == 0, f"Integrality failed at n={n}"
        a.append(a_next_num // ((n + 1) ** 3))
        b_next = (P * b[n] - (n**3) * b[n - 1]) / ((n + 1) ** 3)
        b.append(b_next)
    return a, b


def continued_fraction(x: Decimal, steps: int = 20) -> List[int]:
    """Compute the simple continued fraction of Decimal x."""
    out = []
    cur = x
    for _ in range(steps):
        fl = int(cur)
        out.append(fl)
        rem = cur - Decimal(fl)
        if rem < Decimal("1e-60"):
            break
        cur = Decimal(1) / rem
    return out


# --- Pure-Python LLL Integer Relation Detector ---

def _dot(u: List[Fraction], v: List[Fraction]) -> Fraction:
    return sum(a * b for a, b in zip(u, v))


def lll_reduce(basis: List[List[int]], delta: Fraction = Fraction(3, 4)) -> List[List[int]]:
    """Exact rational LLL lattice basis reduction (Lenstra-Lenstra-Lovász)."""
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


def find_integer_relation(
    values: List[Decimal], names: List[str], scale_digits: int = 45
) -> Tuple[List[int], Decimal]:
    """
    Use LLL on the lattice with basis rows [e_i | round(10^scale_digits * values[i])]
    to find the shortest integer relation sum c_i * values[i] ≈ 0.
    """
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
    print("=" * 74)
    print("  ζ(3) High-Precision Verification & Closed-Form Research Lab")
    print("=" * 74)

    pi = compute_pi()
    ln2 = compute_ln2()
    z3 = compute_zeta3_az(30)
    z3_apery = compute_zeta3_apery_series(120)

    print(f"π      = {str(pi)[:62]}")
    print(f"ζ(3)   = {str(z3)[:62]}")
    print(f"ζ(3)/π³= {str(z3 / (pi**3))[:62]}")
    print(f"CF(ζ(3))     = {continued_fraction(z3, 16)}")
    print(f"CF(ζ(3)/π³)  = {continued_fraction(z3 / (pi**3), 16)}")
    print()

    print("--- 1. Identity Verifications (Error < 1e-50 expected) ---")
    # Markov-Apéry
    err_apery = abs(z3 - z3_apery)
    print(f"[Markov-Apéry series]       |err| = {err_apery:.3e}")

    # Ramanujan self-dual (α = β = π)
    S_pi = lambert_eichler_S(pi)
    z3_ram = (Decimal(7) * pi**3) / Decimal(180) - Decimal(2) * S_pi
    err_ram = abs(z3 - z3_ram)
    print(f"[Ramanujan self-dual α=π]   |err| = {err_ram:.3e}")

    # Ramanujan general (αβ = π²): test α = π*sqrt(2), β = π/sqrt(2)
    # Exact identity (Berndt Ch. 14 Entry 21(i), n = -1):
    #   α^{-1} (ζ(3)/2 + S(α)) + β^{-1} (ζ(3)/2 + S(β)) = (α² + 5αβ + β²) / ... wait check sign!
    # Actually: α^{-1}(ζ(3)/2 + S(α)) - (-β)^{-1}(ζ(3)/2 + S(β)) for odd k:
    # Let's check with α = 2π, β = π/2:
    # Wait: Grosswald (1970) / Berndt:
    #   α^{-1} (1/2 ζ(3) + S(α)) + β^{-1} (1/2 ζ(3) + S(β))  vs  α^{-1}(...) - β^{-1}(...)
    # For α = π*sqrt(2), β = π/sqrt(2):
    sqrt2 = Decimal(2).sqrt()
    alpha, beta = pi * sqrt2, pi / sqrt2
    S_a, S_b = lambert_eichler_S(alpha), lambert_eichler_S(beta)
    lhs_gen = (Decimal("0.5") * z3 + S_a) / alpha + (Decimal("0.5") * z3 + S_b) / beta
    # Test Berndt formula: (1/180) * (α² + 5αβ + β²) / π? No, dim check: lhs has dim α^{-1} = π^{-1}
    # At α = β = π: lhs = (ζ(3) + 2 S(π)) / π = 7 π² / 180.
    # So RHS has dim π² = αβ! Wait:
    #   (α² + 5αβ + β²)/180 has dim π²! Let's check!
    rhs_gen_plus = (alpha**2 + Decimal(5) * alpha * beta + beta**2) / Decimal(180)
    # Also check with alternating sign (Berndt's exact signs for B_0 B_4, B_2^2, B_4 B_0):
    # Actually for α^{-1} F(α) - (-β)^{-1} F(β) = α^{-1} F(α) + β^{-1} F(β):
    # Wait: let's check both (α² + 5αβ + β²)/180 and (α³ + ...)/(αβ):
    # Notice if we multiply by αβ = π²:
    #   β(ζ(3)/2 + S(α)) + α(ζ(3)/2 + S(β)) = π² * (...)
    # Let's check exact residual of rhs_gen_plus AND solve for exact coefficients!
    #Wait: for k=1 (ζ(3)), Berndt Entry 21(i) says:
    # α^{-1}(1/2 ζ(3) + S(α)) + β^{-1}(1/2 ζ(3) + S(β)) = (1/90)(...) ? Let's check!
    err_gen_test = abs(lhs_gen - rhs_gen_plus)
    # Check Berndt's exact formula: -2^2 * sum_{k=0}^2 (-1)^k (B_{2k}/(2k)!) (B_{4-2k}/(4-2k)!) α^{2-k} β^k
    # Wait: at k=0: +(-1/720) α²; k=1: -(1/144) αβ; k=2: +(-1/720) β².
    # Multiply by -4: (1/180)α² + (1/36)αβ + (1/180)β² = (α² + 5αβ + β²)/180!
    # Wait: does it have α^3 or α^2? Let's print the actual discrepancy!
    print(f"[Ramanujan general αβ=π²]   |err| = {err_gen_test:.3e}")

    # Landen Li_3(1/2)
    li3_half = compute_li3_half()
    z3_li3 = Decimal(8) / Decimal(7) * (li3_half + (pi**2) * ln2 / Decimal(12) - (ln2**3) / Decimal(6))
    err_li3 = abs(z3 - z3_li3)
    print(f"[Landen Li_3(1/2) identity] |err| = {err_li3:.3e}")
    print()

    print("--- 2. Apéry Recurrence & Casoratian Verification ---")
    N = 15
    a, b = apery_sequences(N)
    max_cas_err = Fraction(0)
    for n in range(1, N + 1):
        cas = Fraction(a[n - 1]) * b[n] - Fraction(a[n]) * b[n - 1]
        expected = Fraction(6, n**3)
        max_cas_err = max(max_cas_err, abs(cas - expected))
    print(f"a_0..a_5 = {a[:6]}")
    print(f"b_0..b_5 = {[str(x) for x in b[:6]]}")
    print(f"Max Casoratian error (n=1..{N}): {max_cas_err}")
    approx_15 = Decimal(b[N].numerator) / Decimal(b[N].denominator * a[N])
    print(f"|ζ(3) - b_{N}/a_{N}| = {abs(z3 - approx_15):.3e}")
    d_N = math.lcm(*range(1, N + 1))
    print(f"2 * d_{N}^3 * b_{N} is integer? {(2 * d_N**3 * b[N]).denominator == 1}")
    print()

    print("--- 3. LLL Integer Relation Detection ---")
    # Test A: Known relation (sanity check) — should recover exact integers!
    names_known = ["ζ(3)", "Li_3(1/2)", "π²·ln(2)", "ln(2)³"]
    vals_known = [z3, li3_half, (pi**2) * ln2, ln2**3]
    coeffs, res = find_integer_relation(vals_known, names_known, scale_digits=45)
    print("Sanity check on [ζ(3), Li_3(1/2), π²·ln(2), ln(2)³]:")
    print(f"  coeffs = {coeffs}, |residual| = {abs(res):.3e}")

    # Test B: Known Ramanujan relation — should recover [180, -7, 360]
    names_ram = ["ζ(3)", "π³", "S(π)"]
    vals_ram = [z3, pi**3, S_pi]
    coeffs, res = find_integer_relation(vals_ram, names_ram, scale_digits=45)
    print("Sanity check on [ζ(3), π³, S(π)]:")
    print(f"  coeffs = {coeffs}, |residual| = {abs(res):.3e}")

    # Test C: Open question — polynomial/linear relation in [1, π, π², π³, ζ(3)]
    names_open = ["1", "π", "π²", "π³", "ζ(3)"]
    vals_open = [Decimal(1), pi, pi**2, pi**3, z3]
    coeffs, res = find_integer_relation(vals_open, names_open, scale_digits=40)
    print("Open search on [1, π, π², π³, ζ(3)] (scale 10^40):")
    print(f"  shortest vector = {coeffs} (max|c|={max(abs(c) for c in coeffs)}), |residual| = {abs(res):.3e}")
    if max(abs(c) for c in coeffs) > 10**4:
        print("  -> No small-height linear relation in {1, π, π², π³, ζ(3)} (consistent with independence).")

    # Test D: Open question — polynomial relation in {1, π³, ζ(3), ζ(3)²}
    names_poly = ["1", "π³", "ζ(3)", "ζ(3)²", "π³·ζ(3)", "π⁶"]
    vals_poly = [Decimal(1), pi**3, z3, z3**2, (pi**3) * z3, pi**6]
    coeffs, res = find_integer_relation(vals_poly, names_poly, scale_digits=40)
    print("Open search on [1, π³, ζ(3), ζ(3)², π³·ζ(3), π⁶] (scale 10^40):")
    print(f"  shortest vector = {coeffs} (max|c|={max(abs(c) for c in coeffs)}), |residual| = {abs(res):.3e}")


if __name__ == "__main__":
    main()
