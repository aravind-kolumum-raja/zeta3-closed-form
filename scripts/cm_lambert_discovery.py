#!/usr/bin/env python3
"""
cm_lambert_discovery.py — 100-digit LLL discovery engine for exact closed forms
of ζ(3) and π³ in CM quadratic Lambert / Eichler-integral bases.

Tests S_3(r) = sum_{k>=1} 1 / (k^3 * (exp(π r k) - 1)) across imaginary quadratic
CM fields Q(sqrt(-d)), Heegner discriminants, and alternating Lambert series,
discovering exact integer-coefficient closed forms and ranking them by convergence
rate (decimal digits gained per term).
"""

from __future__ import annotations

import math
from decimal import Decimal, getcontext
from fractions import Fraction
from typing import Dict, List, Tuple

WORK_PREC = 110
getcontext().prec = WORK_PREC


def compute_pi() -> Decimal:
    C = 426880 * Decimal(10005).sqrt()
    M, L, X, K, S = Decimal(1), Decimal(13591409), Decimal(1), Decimal(6), Decimal(13591409)
    for i in range(1, 10):
        M = (M * (K**3 - 16 * K)) / Decimal(i**3)
        L += 545140134
        X *= -262537412640768000
        S += (M * L) / X
        K += 12
    return C / S


def compute_zeta3() -> Decimal:
    """Wedeniwski (1998) WZ series: ~5.04 digits/term -> 24 terms gives >115 digits."""
    s = Decimal(0)
    for n in range(24):
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


def S3(r: Decimal, pi: Decimal, sign_den: int = -1, alt_num: bool = False) -> Decimal:
    """
    Generalized weight-3 Lambert/Eichler series:
      sum_{k>=1} ( (-1)^(k-1) if alt_num else 1 ) / (k^3 * (exp(π r k) + sign_den))
    """
    s = Decimal(0)
    k = 1
    while True:
        ex = (pi * r * Decimal(k)).exp()
        num = Decimal((-1) ** (k - 1)) if alt_num else Decimal(1)
        term = num / (Decimal(k**3) * (ex + Decimal(sign_den)))
        s += term
        if abs(term) < Decimal("1e-106"):
            break
        k += 1
    return s


def _dot(u: List[Fraction], v: List[Fraction]) -> Fraction:
    return sum(a * b for a, b in zip(u, v))


def lll_reduce(basis: List[List[int]], delta: Fraction = Fraction(3, 4)) -> List[List[int]]:
    n = len(basis)
    B = [list(map(int, row)) for row in basis]

    def gs():
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

    B_star, mu, norms = gs()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > Fraction(1, 2):
                q = round(mu[k][j])
                B[k] = [B[k][m] - q * B[j][m] for m in range(len(B[k]))]
                B_star, mu, norms = gs()
        if norms[k] >= (delta - mu[k][k - 1] ** 2) * norms[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            B_star, mu, norms = gs()
            k = max(k - 1, 1)
    return B


def find_relations(values: List[Decimal], scale_digits: int = 80) -> List[Tuple[List[int], Decimal]]:
    """Return ALL exact short integer relations found in the LLL-reduced basis."""
    S = Decimal(10) ** scale_digits
    m = len(values)
    basis = []
    for i in range(m):
        row = [0] * (m + 1)
        row[i] = 1
        row[m] = int((values[i] * S).to_integral_value())
        basis.append(row)
    reduced = lll_reduce(basis)
    hits = []
    for row in reduced:
        coeffs = row[:m]
        if all(c == 0 for c in coeffs):
            continue
        res = sum(Decimal(c) * v for c, v in zip(coeffs, values))
        if abs(res) < Decimal("1e-70") and max(abs(c) for c in coeffs) < 10**7:
            # Canonicalize sign
            for c in coeffs:
                if c != 0:
                    if c < 0:
                        coeffs = [-x for x in coeffs]
                    break
            hits.append((coeffs, res))
    return hits


def format_relation(coeffs: List[int], names: List[str]) -> str:
    parts = []
    for c, nm in zip(coeffs, names):
        if c == 0:
            continue
        sign = "+" if c > 0 else "-"
        ac = abs(c)
        term = f"{nm}" if ac == 1 else f"{ac}*{nm}"
        if not parts:
            parts.append(term if c > 0 else f"-{term}")
        else:
            parts.append(f" {sign} {term}")
    return "".join(parts) + " = 0"


def main() -> None:
    print("=" * 80)
    print("  100-Digit CM Quadratic Lambert / Eichler Closed-Form Discovery Engine")
    print("=" * 80)
    pi = compute_pi()
    pi3 = pi**3
    z3 = compute_zeta3()
    print(f"ζ(3) (100d) = {str(z3)[:72]}...")
    print()

    # 1. Systematic sweep over quadratic & rational modular curves:
    #    For each discriminant d, test Ramanujan-Plouffe modular ladders:
    #    r_set = {2/sqrt(d), 2*sqrt(d), 4/sqrt(d), 4*sqrt(d)} and {1/sqrt(d), sqrt(d), ...}
    print("--- 1. Exact Closed Forms for ζ(3) Alone (Eliminating π³!) Across CM Fields ---")
    found_z3_formulas: List[Tuple[float, str, Decimal]] = []

    # Test integer level N families: S_3(2/N), S_3(2), S_3(2N) and alternating partners
    for N in [2, 3, 4, 5, 6, 7]:
        r_vals = [Decimal(2) / Decimal(N), Decimal(2), Decimal(2 * N)]
        names = ["ζ(3)", "π³", f"S3(2/{N})", "S3(2)", f"S3({2*N})"]
        vals = [z3, pi3] + [S3(r, pi) for r in r_vals]
        # First check pure ζ(3) without π³!
        hits_pure = find_relations([z3] + vals[2:], scale_digits=85)
        for c, res in hits_pure:
            if c[0] != 0:
                min_r = float(min(r_vals[i] for i in range(len(r_vals)) if c[i + 1] != 0))
                digits_per_term = min_r * math.pi / math.log(10)
                found_z3_formulas.append((digits_per_term, format_relation(c, ["ζ(3)"] + names[2:]), abs(res)))
        # Also check with π³ included (to discover clean 2-term/3-term modular relations)
        hits_all = find_relations(vals, scale_digits=85)
        for c, res in hits_all:
            if c[0] != 0 and c[1] != 0:
                print(f"  [Level N={N} with π³] {format_relation(c, names)}  (|res|={abs(res):.1e})")

    # Test CM imaginary quadratic fields Q(sqrt(-d)): r = 2/sqrt(d), 2*sqrt(d), 4/sqrt(d), 4*sqrt(d)
    print("\n--- 2. Imaginary Quadratic CM Discriminants Q(sqrt(-d)) (incl. Heegner) ---")
    for d in [2, 3, 5, 6, 7, 10, 11, 13, 15, 19, 43, 67, 163]:
        sd = Decimal(d).sqrt()
        # Self-dual pair under alpha*beta = pi^2 is (alpha, beta) = (pi*sqrt(d), pi/sqrt(d)),
        # corresponding to S_3(2*sqrt(d)) and S_3(2/sqrt(d))!
        s_hi = S3(Decimal(2) * sd, pi)
        s_lo = S3(Decimal(2) / sd, pi)
        # Note: Ramanujan's general formula relates ζ(3), π³*sqrt(d), s_hi, s_lo!
        # Let's test BOTH [ζ(3), π³*sqrt(d), s_hi, s_lo] AND pure [ζ(3), s_hi, s_lo, S3(2)]!
        names_cm = ["ζ(3)", f"π³√{d}", f"S3(2√{d})", f"S3(2/√{d})"]
        vals_cm = [z3, pi3 * sd, s_hi, s_lo]
        hits = find_relations(vals_cm, scale_digits=85)
        for c, res in hits:
            if c[0] != 0:
                print(f"  [CM d={d:3d}] {format_relation(c, names_cm)}  (|res|={abs(res):.1e})")

        # Check eliminating π³ using S3(2) + S3(2√d) + S3(2/√d) when √d is rational? No, π³√d vs π³!
        # Instead eliminate π³√d using TWO pairs at the SAME d: (2√d, 2/√d) and (4√d, 1/√d) -> ALL have π³√d!
        s_hi2 = S3(Decimal(4) * sd, pi)
        s_lo2 = S3(Decimal(1) / sd, pi)
        names_el = ["ζ(3)", f"S3(1/√{d})", f"S3(2/√{d})", f"S3(2√{d})", f"S3(4√{d})"]
        vals_el = [z3, s_lo2, s_lo, s_hi, s_hi2]
        hits_el = find_relations(vals_el, scale_digits=85)
        for c, res in hits_el:
            if c[0] != 0:
                min_r = float(min([Decimal(1)/sd, Decimal(2)/sd, Decimal(2)*sd, Decimal(4)*sd][i] for i in range(4) if c[i+1] != 0))
                dpt = min_r * math.pi / math.log(10)
                found_z3_formulas.append((dpt, format_relation(c, names_el), abs(res)))

    # 3. Test alternating Lambert series A_3^-(r), S_3^+(r), A_3^+(r) at r = 1, 2, 3, 4!
    print("\n--- 3. Alternating & Plus-Sign Lambert Series (Fastest Convergence Hunt!) ---")
    basis_alt_names = [
        "ζ(3)", "π³",
        "S3-(1)", "S3-(2)", "S3-(4)",
        "S3+(1)", "S3+(2)",
        "A3-(1)", "A3-(2)",
        "A3+(1)", "A3+(2)",
    ]
    basis_alt_vals = [
        z3, pi3,
        S3(Decimal(1), pi, -1, False),
        S3(Decimal(2), pi, -1, False),
        S3(Decimal(4), pi, -1, False),
        S3(Decimal(1), pi, +1, False),
        S3(Decimal(2), pi, +1, False),
        S3(Decimal(1), pi, -1, True),
        S3(Decimal(2), pi, -1, True),
        S3(Decimal(1), pi, +1, True),
        S3(Decimal(2), pi, +1, True),
    ]
    # Check small clean subsets for ζ(3) without π³ AND with π³
    for idx_sub, tag in [
        ([0, 1, 6], "ζ(3), π³, S3+(2)"),
        ([0, 1, 8], "ζ(3), π³, A3-(2)"),
        ([0, 1, 10], "ζ(3), π³, A3+(2)"),
        ([0, 1, 5, 6], "ζ(3), π³, S3+(1), S3+(2)"),
        ([0, 1, 7, 8], "ζ(3), π³, A3-(1), A3-(2)"),
        ([0, 1, 9, 10], "ζ(3), π³, A3+(1), A3+(2)"),
        ([0, 3, 6, 8, 10], "Pure r=2: ζ(3), S3-(2), S3+(2), A3-(2), A3+(2)"),
        ([0, 1, 3, 6, 8, 10], "r=2 + π³"),
    ]:
        sub_n = [basis_alt_names[i] for i in idx_sub]
        sub_v = [basis_alt_vals[i] for i in idx_sub]
        hits = find_relations(sub_v, scale_digits=85)
        for c, res in hits:
            if c[0] != 0:
                print(f"  [{tag}] {format_relation(c, sub_n)}  (|res|={abs(res):.1e})")

    # 4. Print ranked pure-ζ(3) Lambert formulas (no π³!)
    print("\n--- 4. Ranked Pure-ζ(3) Closed Forms in Lambert/Eichler Bases (No π³!) ---")
    found_z3_formulas.sort(key=lambda x: -x[0])
    seen = set()
    for dpt, expr, res in found_z3_formulas:
        if expr in seen:
            continue
        seen.add(expr)
        print(f"  [{dpt:.3f} digits/term] {expr}  (|res|={res:.1e})")

    # 5. High-precision (90-digit) exclusion check on ζ(3)/π³ itself!
    print("\n--- 5. 95-Digit Continued Fraction & Height Exclusion for ζ(3)/π³ ---")
    ratio = z3 / pi3
    cf = []
    cur = ratio
    for _ in range(50):
        fl = int(cur)
        cf.append(fl)
        rem = cur - Decimal(fl)
        if rem < Decimal("1e-90"):
            break
        cur = Decimal(1) / rem
    print(f"  First 50 CF partial quotients of ζ(3)/π³:\n    {cf}")
    print(f"  Max partial quotient in first 50: {max(cf)} (generic Gauss-Kuzmin behavior)")
    # Compute 50th convergent denominator (certifies no exact rational p/q with q < q_50!)
    p_prev, p_cur = 0, 1
    q_prev, q_cur = 1, 0
    for a_i in cf[:45]:
        p_prev, p_cur = p_cur, a_i * p_cur + p_prev
        q_prev, q_cur = q_cur, a_i * q_cur + q_prev
    print(f"  45th convergent denominator q_45 = {q_cur} (~10^{math.log10(q_cur):.1f})")
    print(f"  => Certifies ζ(3)/π³ ≠ p/q for ALL integers |q| < 10^{math.log10(q_cur):.1f}!")


if __name__ == "__main__":
    main()
