#!/usr/bin/env python3
"""
225-Digit Independent Certifier for All 8 Pure (Zero-pi, Zero-ln)
Algebraic-Unit Trilogarithm Closed Forms for zeta(3)
=================================================================
Verifies to 225 decimal digits all 8 pure positive-power Li_3 closed forms
discovered in our overnight campaign across:
  - All 4 trinomial Pisot numbers x^p - x^{p-1} - 1 = 0 (p = 2, 3, 4, 5):
      1. Golden ratio phi              (x^2 - x - 1 = 0)
      2. Supergolden psi / Weber D=-31 (x^3 - x^2 - 1 = 0)
      3. 2nd smallest Pisot theta_2    (x^4 - x^3 - 1 = 0)
      4. Plastic rho / Weber D=-23     (x^3 - x - 1 = 0 | x^5 - x^4 - 1)
  - All 3 smallest members (a = 1, 2, 3) of our infinite Salem collapse family
    (1 - x^{-a})(1 - x^{-(a+2)}) = x^{-(a+1)}:
      5. Smallest degree-4 Salem sigma_4 (x^4 - x^3 - x^2 - x + 1 = 0) -- 1*zeta(3)!
      6. Smallest degree-6 Salem sigma_6 (x^6 - x^4 - x^3 - x^2 + 1 = 0)
      7. Smallest degree-8 Salem sigma_8 (x^8 - x^5 - x^4 - x^3 + 1 = 0)
  - Lehmer's degree-10 Salem number (smallest known Mahler measure > 1!):
      8. Lehmer Salem lambda_L (x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1 = 0)
"""

from decimal import Decimal, getcontext

getcontext().prec = 235


def compute_zeta3():
    s = Decimal(0)
    b = 1
    for k in range(1, 410):
        b = b * (2 * k - 1) * (2 * k) // (k * k)
        term = Decimal(1) / Decimal(k ** 3 * b)
        s += term if (k % 2 == 1) else -term
    return s * Decimal(5) / Decimal(2)


def solve_root(coeffs, x0):
    x = Decimal(str(x0))
    d = len(coeffs) - 1
    for _ in range(55):
        v = Decimal(0)
        dv = Decimal(0)
        for i, c in enumerate(coeffs):
            p = d - i
            v = v * x + c
            if p > 0:
                dv = dv * x + c * p
        dx = v / dv
        x -= dx
        if abs(dx) < Decimal(10) ** (-230):
            break
    # Verify polynomial residual
    v = Decimal(0)
    for c in coeffs:
        v = v * x + c
    assert abs(v) < Decimal(10) ** (-225)
    return x


def li3(x):
    eps = Decimal(10) ** (-231)
    s = Decimal(0)
    xp = x
    n = 1
    while True:
        t = xp / Decimal(n * n * n)
        s += t
        if abs(t) < eps:
            break
        xp *= x
        n += 1
    return s


PURE_LADDERS = [
    ("1. Golden ratio phi (x^2-x-1)", [1, -1, -1], 1.618033989, 3,
     [(1, 12), (2, -6), (3, -16), (4, 9), (6, 8), (12, -1)]),
    ("2. Supergolden psi / Weber D=-31 (x^3-x^2-1)", [1, -1, 0, -1], 1.465571232, 12,
     [(1, 52), (2, -58), (3, 8), (4, -10), (5, 16), (6, 11), (10, -4)]),
    ("3. 2nd smallest Pisot theta_2 (x^4-x^3-1)", [1, -1, 0, 0, -1], 1.380277569, 12,
     [(1, 24), (2, 31), (3, -64), (4, 2), (5, 8), (6, -20), (7, 28), (9, 16), (14, -7)]),
    ("4. Plastic rho / Weber D=-23 (x^3-x-1)", [1, 0, -1, -1], 1.324717957, 4,
     [(1, 12), (2, -5), (3, -4), (4, -8), (5, 8), (8, 2)]),
    ("5. Smallest deg-4 Salem sigma_4 (1*zeta(3)!)", [1, -1, -1, -1, 1], 1.722083806, 1,
     [(1, 5), (2, -7), (3, 2), (4, 4), (5, -8), (6, 2), (8, 8), (10, 1), (11, -1), (12, -2)]),
    ("6. Smallest deg-6 Salem sigma_6", [1, 0, -1, -1, -1, 0, 1], 1.401268368, 21,
     [(1, 58), (2, 15), (3, -106), (4, 35), (5, -24), (6, 14), (7, 56), (9, -14), (10, 3), (11, -10), (12, 1)]),
    ("7. Smallest deg-8 Salem sigma_8", [1, 0, 0, -1, -1, -1, 0, 0, 1], 1.280638156, 40,
     [(1, 112), (2, -57), (3, 2), (4, -46), (5, 14), (6, -65), (7, 66), (8, -36), (9, 64), (10, 36),
      (12, -5), (13, -18), (14, -41), (16, 9), (18, 9)]),
    ("8. Lehmer deg-10 Salem lambda_L (min M>1!)", [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], 1.176280818, 108,
     [(1, 194), (2, -99), (3, 224), (4, -211), (5, 126), (6, -224), (7, 116), (8, -183), (9, 4),
      (10, -94), (11, 186), (12, 59), (13, 114), (16, -34), (17, -16), (18, -57), (20, 17)]),
]


def main():
    Z3 = compute_zeta3()
    print("=" * 78)
    print("225-DIGIT CERTIFICATION OF ALL 8 PURE ZERO-PI ZERO-LN Li_3 FORMULAS FOR zeta(3)")
    print("=" * 78)
    for name, poly, x0, c0, terms in PURE_LADDERS:
        a = solve_root(poly, x0)
        inv = Decimal(1) / a
        rhs = sum(Decimal(c) * li3(inv ** k) for k, c in terms)
        err = abs(Decimal(c0) * Z3 - rhs)
        assert err < Decimal(10) ** (-220), f"Failed {name}: {err}"
        print(f"  [PASS {float(err):.2e}] {name:46s} | a={float(a):.10f} | {len(terms):2d} terms, c0={c0}")
    print("\nALL 8 PURE CLOSED FORMS VERIFIED TO > 220 DECIMAL DIGITS!")


if __name__ == "__main__":
    main()
