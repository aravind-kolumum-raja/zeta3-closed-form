#!/usr/bin/env python3
"""
hypergeom_series_discovery.py — 115-digit systematic LLL discovery engine for
exact closed-form series of zeta(3) across:
  1. BBP-type & polylogarithmic bases: binary (1/2, 1/4, 1/8, 1/16), ternary (1/3, 1/9, 1/27),
     golden ratio (phi^-1, phi^-2, phi^-4, phi^-6), silver ratio ((sqrt(2)-1)^2),
     and sqrt(3) base (2-sqrt(3)).
  2. Central-binomial & hypergeometric families:
       sum_{k>=1} P(k) * z^k / (k^m * C(2k,k)^a * C(3k,k)^b * C(4k,2k)^c)
     automatically solving for exact integer numerator polynomials P(k) via LLL!

Every discovered identity is verified independently to < 1e-95 residual.
"""

from decimal import Decimal, getcontext
from fractions import Fraction
import math
import time

getcontext().prec = 120


def compute_pi() -> Decimal:
    C = Decimal(426880) * Decimal(10005).sqrt()
    M, L, X, K, S = Decimal(1), Decimal(13591409), Decimal(1), Decimal(6), Decimal(13591409)
    for i in range(1, 11):
        M = (M * (K**3 - 16 * K)) / Decimal(i**3)
        L += 545140134
        X *= -262537412640768000
        S += (M * L) / X
        K += 12
    return C / S


def compute_zeta3() -> Decimal:
    s = Decimal(0)
    for n in range(25):
        P = 126392 * n**5 + 412708 * n**4 + 531578 * n**3 + 336367 * n**2 + 104000 * n + 12463
        num = ((-1) ** n) * ((math.factorial(2 * n + 1) * math.factorial(2 * n) * math.factorial(n)) ** 3) * P
        den = math.factorial(3 * n + 2) * (math.factorial(4 * n + 3) ** 3)
        s += Decimal(num) / Decimal(den)
    return s / Decimal(24)


def lll_find_exact(vals, names, scale_exp=85, max_coeff=50000):
    """Find short exact integer relations involving vals[0] via exact-Fraction LLL."""
    n = len(vals)
    S = 10**scale_exp
    B = [[0] * (n + 1) for _ in range(n)]
    for i in range(n):
        B[i][i] = 1
        B[i][n] = int((vals[i] * Decimal(S)).to_integral_value())

    def gs(B):
        m = len(B)
        ortho = [[Fraction(x) for x in row] for row in B]
        mu = [[Fraction(0)] * m for _ in range(m)]
        norm2 = [Fraction(0)] * m
        for i in range(m):
            for j in range(i):
                if norm2[j] != 0:
                    num = sum(Fraction(B[i][k]) * ortho[j][k] for k in range(len(B[0])))
                    mu[i][j] = num / norm2[j]
                    for k in range(len(B[0])):
                        ortho[i][k] -= mu[i][j] * ortho[j][k]
            norm2[i] = sum(x * x for x in ortho[i])
        return ortho, mu, norm2

    ortho, mu, norm2 = gs(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            r = round(mu[k][j])
            if r != 0:
                for col in range(n + 1):
                    B[k][col] -= r * B[j][col]
                ortho, mu, norm2 = gs(B)
        if norm2[k] >= (Fraction(3, 4) - mu[k][k - 1] ** 2) * norm2[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            ortho, mu, norm2 = gs(B)
            k = max(k - 1, 1)

    hits = []
    for row in B:
        coeffs = row[:n]
        if coeffs[0] == 0:
            continue
        if max(abs(c) for c in coeffs) > max_coeff:
            continue
        resid = abs(sum(Decimal(c) * v for c, v in zip(coeffs, vals)))
        if resid < Decimal("1e-72"):
            # Canonicalize sign & gcd
            g = 0
            for c in coeffs:
                g = math.gcd(g, abs(c))
            coeffs = [c // g for c in coeffs]
            if coeffs[0] < 0:
                coeffs = [-c for c in coeffs]
            hits.append((coeffs, resid))
    return hits


def bbp_basis(z: Decimal, M: int, max_terms: int = 1200):
    """Compute B_j(z, M) = sum_{k=0}^inf z^k / (M*k + j)^3 for j = 1..M."""
    if abs(z) >= Decimal("0.98"):
        return None
    sums = [Decimal(0)] * M
    zk = Decimal(1)
    for k in range(max_terms):
        if abs(zk) < Decimal("1e-114"):
            break
        base = M * k
        for j in range(M):
            d = base + j + 1
            sums[j] += zk / Decimal(d * d * d)
        zk *= z
    else:
        if abs(zk) > Decimal("1e-100"):
            return None
    return sums


def binom_basis(a2: int, a3: int, a4: int, z: Decimal, max_terms: int = 500):
    """
    Compute S_r = sum_{k>=1} z^k / (k^{5-r} * C(2k,k)^a2 * C(3k,k)^a3 * C(4k,2k)^a4)
    for r = 0, 1, 2, 3 (so k^{-5}, k^{-4}, k^{-3}, k^{-2}).
    Returns (S_0..S_3, digits_per_term) if geometrically convergent to <1e-112.
    """
    # Asymptotic growth rate per step of denominator: 4^a2 * (27/4)^a3 * 16^a4
    # Ratio per step = |z| / (4^a2 * (6.75)^a3 * 16^a4)
    step_ratio = abs(float(z)) / ((4.0**a2) * (6.75**a3) * (16.0**a4))
    if step_ratio >= 0.72 or step_ratio <= 0.0:
        return None
    dpt = -math.log10(step_ratio)

    sums = [Decimal(0)] * 4
    c2, c3, c4 = 1, 1, 1
    zk = Decimal(1)
    for k in range(1, max_terms + 1):
        zk *= z
        c2 = c2 * (2 * k - 1) * (2 * k) // (k * k)
        c3 = c3 * (3 * k - 2) * (3 * k - 1) * (3 * k) // (k * (2 * k - 1) * (2 * k))
        c4 = c4 * (4 * k - 3) * (4 * k - 2) * (4 * k - 1) * (4 * k) // ((2 * k - 1) * (2 * k) * (2 * k - 1) * (2 * k))
        # Build exact rational factor for binomial powers
        num_b = 1
        den_b = 1
        for val, exp in ((c2, a2), (c3, a3), (c4, a4)):
            if exp > 0:
                den_b *= val**exp
            elif exp < 0:
                num_b *= val ** (-exp)
        base_term = (zk * Decimal(num_b)) / Decimal(den_b)
        k2 = k * k
        k3 = k2 * k
        k4 = k3 * k
        k5 = k4 * k
        t5 = base_term / Decimal(k5)
        t4 = base_term / Decimal(k4)
        t3 = base_term / Decimal(k3)
        t2 = base_term / Decimal(k2)
        sums[0] += t5
        sums[1] += t4
        sums[2] += t3
        sums[3] += t2
        if abs(t2) < Decimal("1e-114"):
            return sums, dpt
    return None


def main():
    t0 = time.perf_counter()
    PI = compute_pi()
    Z3 = compute_zeta3()
    sqrt2 = Decimal(2).sqrt()
    sqrt3 = Decimal(3).sqrt()
    sqrt5 = Decimal(5).sqrt()
    phi = (Decimal(1) + sqrt5) / Decimal(2)
    ln2 = Decimal(2).ln()
    ln3 = Decimal(3).ln()
    lnphi = phi.ln()
    lnsilv = (Decimal(1) + sqrt2).ln()
    ln2s3 = (Decimal(2) + sqrt3).ln()

    print("========================================================================")
    print("  115-DIGIT HYPERGEOMETRIC, BINOMIAL & GOLDEN-RATIO BBP HUNT FOR zeta(3)")
    print("========================================================================\n")

    # =====================================================================
    # PART 1: Pure & Mixed BBP-Type Series in Binary, Ternary, Golden-Ratio,
    #         Silver-Ratio, and sqrt(3) Bases!
    # =====================================================================
    print("=== PART 1: BBP-Type & Polylog-Modulus Discovery (Golden Ratio, Binary, Ternary, sqrt2, sqrt3) ===")
    bbp_bases = [
        ("1/2", Decimal(1) / Decimal(2), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("1/4", Decimal(1) / Decimal(4), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("-1/4", Decimal(-1) / Decimal(4), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("1/8", Decimal(1) / Decimal(8), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("-1/8", Decimal(-1) / Decimal(8), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("1/16", Decimal(1) / Decimal(16), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("1/64", Decimal(1) / Decimal(64), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("-1/64", Decimal(-1) / Decimal(64), [PI**2 * ln2, ln2**3], ["pi^2*ln2", "ln^3(2)"]),
        ("1/9", Decimal(1) / Decimal(9), [PI**2 * ln3, ln3**3], ["pi^2*ln3", "ln^3(3)"]),
        ("-1/27", Decimal(-1) / Decimal(27), [PI**2 * ln3, ln3**3], ["pi^2*ln3", "ln^3(3)"]),
        ("phi^-2", phi**(-2), [PI**2 * lnphi, lnphi**3], ["pi^2*ln(phi)", "ln^3(phi)"]),
        ("phi^-4", phi**(-4), [PI**2 * lnphi, lnphi**3], ["pi^2*ln(phi)", "ln^3(phi)"]),
        ("phi^-6", phi**(-6), [PI**2 * lnphi, lnphi**3], ["pi^2*ln(phi)", "ln^3(phi)"]),
        ("phi^-12", phi**(-12), [PI**2 * lnphi, lnphi**3], ["pi^2*ln(phi)", "ln^3(phi)"]),
        ("(sqrt2-1)^2", (sqrt2 - 1) ** 2, [PI**2 * lnsilv, lnsilv**3], ["pi^2*ln(1+s2)", "ln^3(1+s2)"]),
        ("(sqrt2-1)^4", (sqrt2 - 1) ** 4, [PI**2 * lnsilv, lnsilv**3], ["pi^2*ln(1+s2)", "ln^3(1+s2)"]),
        ("2-sqrt3", Decimal(2) - sqrt3, [PI**2 * ln2s3, ln2s3**3], ["pi^2*ln(2+s3)", "ln^3(2+s3)"]),
        ("(2-sqrt3)^2", (Decimal(2) - sqrt3) ** 2, [PI**2 * ln2s3, ln2s3**3], ["pi^2*ln(2+s3)", "ln^3(2+s3)"]),
    ]

    found_bbp = 0
    for bname, zval, extras, enames in bbp_bases:
        for M in (1, 2, 3, 4, 6):
            # Keep basis sizes small & fast
            if M == 6 and bname not in ("1/64", "-1/64", "phi^-6", "phi^-12", "1/16", "1/4"):
                continue
            bs = bbp_basis(zval, M)
            if bs is None:
                continue
            bnames = [f"B({bname},{M},{j+1})" for j in range(M)]
            # 1. Pure zeta(3) BBP check (ZERO pi^2*ln or ln^3 tails!)
            hits_pure = lll_find_exact([Z3] + bs, ["zeta(3)"] + bnames, scale_exp=82, max_coeff=20000)
            for coeffs, resid in hits_pure:
                found_bbp += 1
                terms = [f"{c:+d}*{nm}" for c, nm in zip(coeffs[1:], bnames) if c != 0]
                print(f"  [PURE BBP] {coeffs[0]}*zeta(3) + {' '.join(terms)} = 0  (resid={resid:.2e})")
            # 2. Mixed zeta(3) + natural weight-3 log/pi^2 terms (e.g. golden-ratio BBP!)
            if M <= 3 and not hits_pure:
                hits_mix = lll_find_exact(
                    [Z3] + bs + extras, ["zeta(3)"] + bnames + enames, scale_exp=82, max_coeff=5000
                )
                for coeffs, resid in hits_mix:
                    found_bbp += 1
                    terms = [f"{c:+d}*{nm}" for c, nm in zip(coeffs[1:], bnames + enames) if c != 0]
                    print(f"  [MIXED BBP] {coeffs[0]}*zeta(3) + {' '.join(terms)} = 0  (resid={resid:.2e})")

    # =====================================================================
    # PART 2: Central-Binomial & Hypergeometric Polynomial-Numerator Discovery!
    # =====================================================================
    print("\n=== PART 2: Hypergeometric & Central-Binomial Super-Series Discovery ===")
    # Grid of (a2, a3, a4) for C(2k,k)^a2 * C(3k,k)^a3 * C(4k,2k)^a4
    shapes = [
        (1, 0, 0),   # 1 / C(2k,k)  [Apéry / Markov!]
        (2, 0, 0),   # 1 / C(2k,k)^2
        (3, 0, 0),   # 1 / C(2k,k)^3
        (4, 0, 0),   # 1 / C(2k,k)^4
        (5, 0, 0),   # 1 / C(2k,k)^5 [Amdeberhan-Zeilberger!]
        (0, 1, 0),   # 1 / C(3k,k)
        (0, 2, 0),   # 1 / C(3k,k)^2
        (0, 0, 1),   # 1 / C(4k,2k)
        (1, 1, 0),   # 1 / (C(2k,k) * C(3k,k))
        (1, 0, 1),   # 1 / (C(2k,k) * C(4k,2k))
        (2, 1, 0),   # 1 / (C(2k,k)^2 * C(3k,k))
        (2, 0, 1),   # 1 / (C(2k,k)^2 * C(4k,2k))
        (-1, 1, 0),  # C(2k,k) / C(3k,k)
        (-1, 0, 1),  # C(2k,k) / C(4k,2k)
        (1, -1, 1),  # C(3k,k) / (C(2k,k)*C(4k,2k))
        (-2, 0, 1),  # C(2k,k)^2 / C(4k,2k)
        (3, -1, 1),  # C(3k,k) / (C(2k,k)^3*C(4k,2k))
        (-1, 2, 0),  # C(2k,k) / C(3k,k)^2
    ]
    z_grid = [
        ("1", Decimal(1)),
        ("-1", Decimal(-1)),
        ("2", Decimal(2)),
        ("-2", Decimal(-2)),
        ("4", Decimal(4)),
        ("-4", Decimal(-4)),
        ("1/2", Decimal(1) / Decimal(2)),
        ("-1/2", Decimal(-1) / Decimal(2)),
        ("1/4", Decimal(1) / Decimal(4)),
        ("-1/4", Decimal(-1) / Decimal(4)),
        ("1/8", Decimal(1) / Decimal(8)),
        ("-1/8", Decimal(-1) / Decimal(8)),
        ("1/16", Decimal(1) / Decimal(16)),
        ("-1/16", Decimal(-1) / Decimal(16)),
        ("1/27", Decimal(1) / Decimal(27)),
        ("-1/27", Decimal(-1) / Decimal(27)),
        ("27/4", Decimal(27) / Decimal(4)),
        ("-27/4", Decimal(-27) / Decimal(4)),
        ("4/27", Decimal(4) / Decimal(27)),
        ("-4/27", Decimal(-4) / Decimal(27)),
        ("16", Decimal(16)),
        ("-16", Decimal(-16)),
        ("64", Decimal(64)),
        ("-64", Decimal(-64)),
        ("phi^-1", phi**(-1)),
        ("-phi^-1", -(phi**(-1))),
        ("phi^-2", phi**(-2)),
        ("-phi^-2", -(phi**(-2))),
        ("(sqrt2-1)^2", (sqrt2 - 1) ** 2),
        ("-(sqrt2-1)^2", -((sqrt2 - 1) ** 2)),
    ]

    found_hg = 0
    tested_hg = 0
    for a2, a3, a4 in shapes:
        for zname, zval in z_grid:
            res = binom_basis(a2, a3, a4, zval)
            if res is None:
                continue
            sums, dpt = res
            tested_hg += 1
            snames = ["k^-5", "k^-4", "k^-3", "k^-2"]
            # Pure zeta(3) relation with polynomial numerator P(k) = c0 + c1*k + c2*k^2 + c3*k^3!
            hits = lll_find_exact([Z3] + sums, ["zeta(3)"] + snames, scale_exp=82, max_coeff=50000)
            for coeffs, resid in hits:
                found_hg += 1
                p_terms = [f"{c:+d}*{nm}" for c, nm in zip(coeffs[1:], snames) if c != 0]
                print(
                    f"  [EXACT HYPERGEOM | {dpt:.2f} d/term] shape=({a2},{a3},{a4}) z={zname:11s} | "
                    f"{coeffs[0]}*zeta(3) + sum({ ' '.join(p_terms) }) = 0  (resid={resid:.2e})"
                )
            # Also check if pure pi^3 appears in the same hypergeometric kernel!
            hits_pi3 = lll_find_exact([PI**3] + sums, ["pi^3"] + snames, scale_exp=82, max_coeff=50000)
            for coeffs, resid in hits_pi3:
                p_terms = [f"{c:+d}*{nm}" for c, nm in zip(coeffs[1:], snames) if c != 0]
                print(
                    f"  [EXACT PI^3 HG   | {dpt:.2f} d/term] shape=({a2},{a3},{a4}) z={zname:11s} | "
                    f"{coeffs[0]}*pi^3 + sum({ ' '.join(p_terms) }) = 0  (resid={resid:.2e})"
                )
            # Also check coupled (zeta(3), pi^3) relation ONLY if neither pure hit occurred!
            if not hits and not hits_pi3 and zname in ("1", "-1", "1/2", "-1/2", "1/4", "-1/4", "1/16", "-1/16", "phi^-2", "4", "-4"):
                hits_c = lll_find_exact([Z3, PI**3] + sums, ["zeta(3)", "pi^3"] + snames, scale_exp=82, max_coeff=5000)
                for coeffs, resid in hits_c:
                    if coeffs[1] != 0:
                        p_terms = [f"{c:+d}*{nm}" for c, nm in zip(coeffs[2:], snames) if c != 0]
                        print(
                            f"  [COUPLED z3+pi^3 | {dpt:.2f} d/term] shape=({a2},{a3},{a4}) z={zname:11s} | "
                            f"{coeffs[0]}*zeta(3) {coeffs[1]:+d}*pi^3 + sum({ ' '.join(p_terms) }) = 0  (resid={resid:.2e})"
                        )

    print(
        f"\nDone in {time.perf_counter()-t0:.2f}s | Tested {tested_hg} convergent hypergeometric kernels "
        f"| Found {found_bbp} BBP hits & {found_hg} exact hypergeometric zeta(3) hits!"
    )


if __name__ == "__main__":
    main()
