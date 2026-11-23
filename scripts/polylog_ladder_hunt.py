#!/usr/bin/env python3
"""
polylog_ladder_hunt.py — 115-digit LLL discovery of exact pure and mixed
trilogarithm ladders (Lewin-Zagier / Broadhurst BBP-type closed forms) for zeta(3)
in algebraic unit bases:
  - Golden ratio phi = (1 + sqrt(5)) / 2       (phi^2 = phi + 1)
  - Silver ratio 1 + sqrt(2)                   (u^2 = 2u + 1)
  - sqrt(3) unit 2 + sqrt(3)                   (u^2 = 4u - 1)
  - Plastic constant rho ~ 1.324717957...      (rho^3 = rho + 1)
  - Binary base 2 & Ternary base 3 (Broadhurst-type pure BBP checks)
"""

from decimal import Decimal, getcontext
from fractions import Fraction
import math

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


def li3(x: Decimal) -> Decimal:
    s = Decimal(0)
    xp = x
    for n in range(1, 1500):
        t = xp / Decimal(n**3)
        s += t
        if abs(t) < Decimal("1e-116"):
            break
        xp *= x
    return s


def lll_all(vals, scale_exp=85, max_c=50000):
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
                    mu[i][j] = sum(Fraction(B[i][k]) * ortho[j][k] for k in range(len(B[0]))) / norm2[j]
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

    out = []
    for row in B:
        cs = row[:n]
        if cs[0] == 0 or max(abs(c) for c in cs) > max_c:
            continue
        resid = abs(sum(Decimal(c) * v for c, v in zip(cs, vals)))
        if resid < Decimal("1e-75"):
            g = 0
            for c in cs:
                g = math.gcd(g, abs(c))
            cs = [c // g for c in cs]
            if cs[0] < 0:
                cs = [-c for c in cs]
            out.append((cs, resid))
    return out


def main():
    PI = compute_pi()
    Z3 = compute_zeta3()
    phi = (Decimal(1) + Decimal(5).sqrt()) / Decimal(2)
    silv = Decimal(1) + Decimal(2).sqrt()
    r3 = Decimal(2) + Decimal(3).sqrt()
    # Plastic constant rho: real root of rho^3 - rho - 1 = 0
    rho = Decimal("1.32471795724474602596")
    for _ in range(12):
        rho -= (rho**3 - rho - 1) / (Decimal(3) * rho * rho - 1)
    # Supergolden ratio psi: real root of psi^3 - psi^2 - 1 = 0
    sphi = Decimal("1.46557123187676802665")
    for _ in range(12):
        sphi -= (sphi**3 - sphi**2 - 1) / (Decimal(3) * sphi * sphi - Decimal(2) * sphi)

    ladders = [
        ("GOLDEN RATIO phi = (1+sqrt5)/2 [phi^2 = phi+1]", phi, [1, 2, 3, 4, 6, 8, 10, 12]),
        ("SILVER RATIO 1+sqrt2 [u^2 = 2u+1]", silv, [1, 2, 3, 4, 6, 8]),
        ("SQRT(3) UNIT 2+sqrt3 [u^2 = 4u-1]", r3, [1, 2, 3, 4, 6]),
        ("PLASTIC CONSTANT rho [rho^3 = rho+1]", rho, [1, 2, 3, 4, 5, 6, 8, 9, 10, 12]),
        ("SUPERGOLDEN RATIO psi [psi^3 = psi^2+1]", sphi, [1, 2, 3, 4, 5, 6, 8, 9, 12]),
        ("BINARY BASE 2", Decimal(2), [1, 2, 3, 4, 6]),
        ("TERNARY BASE 3", Decimal(3), [1, 2, 3, 4]),
    ]

    print("========================================================================")
    print("  115-DIGIT TRILOGARITHM LADDER HUNT FOR zeta(3) IN ALGEBRAIC BASES")
    print("========================================================================")

    for title, u, exps in ladders:
        print(f"\n=== {title} ===")
        lnu = u.ln()
        vals = [Z3]
        names = ["zeta(3)"]
        for e in exps:
            vals.append(li3(u**(-e)))
            names.append(f"Li3(u^-{e})")
            if e <= 3:
                vals.append(li3(-(u**(-e))))
                names.append(f"Li3(-u^-{e})")
        # 1. Pure zeta(3) ladder (ZERO pi^2*ln(u), ZERO ln^3(u)!)
        pure = lll_all(vals, scale_exp=85, max_c=50000)
        if pure:
            for cs, resid in pure[:3]:
                terms = " ".join(f"{c:+d}*{nm}" for c, nm in zip(cs[1:], names[1:]) if c != 0)
                print(f"  [PURE LADDER!!] {cs[0]}*zeta(3) {terms} = 0  (resid={resid:.2e})")
        else:
            print("  No pure zeta(3) ladder (without pi^2*ln(u) / ln^3(u))")
        # 2. Mixed ladder including pi^2*ln(u) and ln^3(u)
        ext_vals = vals + [PI**2 * lnu, lnu**3]
        ext_names = names + ["pi^2*ln(u)", "ln^3(u)"]
        mix = lll_all(ext_vals, scale_exp=85, max_c=50000)
        for cs, resid in mix[:3]:
            terms = " ".join(f"{c:+d}*{nm}" for c, nm in zip(cs[1:], ext_names[1:]) if c != 0)
            print(f"  [LADDER]        {cs[0]}*zeta(3) {terms} = 0  (resid={resid:.2e})")

    # Verify the simplified positive-power-only forms (eliminating Li3(-x) = 1/4*Li3(x^2) - Li3(x))
    print("\n=== SIMPLIFIED POSITIVE-POWER CLOSED FORMS (115-DIGIT CHECK) ===")
    g_rhs = (
        12 * li3(phi**(-1))
        - 6 * li3(phi**(-2))
        - 16 * li3(phi**(-3))
        + 9 * li3(phi**(-4))
        + 8 * li3(phi**(-6))
        - li3(phi**(-12))
    ) / Decimal(3)
    err_g = abs(Z3 - g_rhs)
    assert err_g < Decimal("1e-114"), f"Golden pure ladder failed: {err_g}"
    print(f"  [PASS] Pure Golden Ratio (zero pi, zero ln):    err = {err_g:.2e}")

    p_rhs = (
        12 * li3(rho**(-1))
        - 5 * li3(rho**(-2))
        - 4 * li3(rho**(-3))
        - 8 * li3(rho**(-4))
        + 8 * li3(rho**(-5))
        + 2 * li3(rho**(-8))
    ) / Decimal(4)
    err_p = abs(Z3 - p_rhs)
    assert err_p < Decimal("1e-114"), f"Plastic pure ladder failed: {err_p}"
    print(f"  [PASS] Pure Plastic Constant (zero pi, zero ln): err = {err_p:.2e}")

    ln3 = Decimal(3).ln()
    t_rhs = (36 * li3(Decimal(1) / Decimal(3)) - 3 * li3(Decimal(1) / Decimal(9)) + 2 * PI**2 * ln3 - 2 * ln3**3) / Decimal(26)
    err_t = abs(Z3 - t_rhs)
    assert err_t < Decimal("1e-114"), f"Ternary ladder failed: {err_t}"
    print(f"  [PASS] Ternary Base-3 (1/3, 1/9 + pi^2 ln3):    err = {err_t:.2e}")


if __name__ == "__main__":
    main()
