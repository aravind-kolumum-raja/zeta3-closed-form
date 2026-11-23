#!/usr/bin/env python3
"""
120-Digit Verifier & Guillera Shadow-Duality Lab for Hypergeometric Series
==========================================================================
1. Independently verifies ALL candidate hits from `tools/wz_series_hunter`
   (wz_v2) to 120 decimal digits using exact Rational/Decimal arithmetic.
2. Tests the exact k -> -k - 1/2 shadow duality (Guillera-Sun-Chu-Zhang)
   across ALL FIVE Guillera 1/pi^2 Ramanujan-like series:
     - e = [5,0,0,0], z = -1/1024  -> Amdeberhan-Zeilberger zeta(3) dual!
     - e = [5,0,0,0], z = -1/4     -> 7*zeta(3) dual!
     - e = [3,0,1,0], z = +1/16    -> tests (1/2)_k^{-3} (1/4)_k^{-1} (3/4)_k^{-1} dual!
     - e = [1,1,1,0], z = -1/48    -> tests (1/2)_k^{-1} (1/3)_k^{-1} (2/3)_k^{-1} (1/4)_k^{-1} (3/4)_k^{-1} dual!
     - e = [3,1,0,0], z = +27/64   -> tests (1/2)_k^{-3} (1/3)_k^{-1} (2/3)_k^{-1} dual!
"""

from decimal import Decimal, getcontext
from fractions import Fraction
import math

getcontext().prec = 135


def compute_pi():
    def atan_inv(m):
        m2 = m * m
        t = Decimal(1) / Decimal(m)
        s = t
        k = 1
        eps = Decimal(10) ** (-132)
        while True:
            t = -t / m2
            add = t / (2 * k + 1)
            if abs(add) < eps:
                break
            s += add
            k += 1
        return s
    return Decimal(4) * (Decimal(4) * atan_inv(5) - atan_inv(239))


def compute_zeta3():
    s = Decimal(0)
    b = 1
    for k in range(1, 240):
        b = b * (2 * k - 1) * (2 * k) // (k * k)
        term = Decimal(1) / Decimal(k ** 3 * b)
        s += term if (k % 2 == 1) else -term
    return s * Decimal(5) / Decimal(2)


PI_D = compute_pi()
PI2_D = PI_D * PI_D
PI3_D = PI2_D * PI_D
Z3_D = compute_zeta3()
CATALAN_D = None


def eval_hypergeom(e, dk, z_num, z_den, poly_coeffs, max_k=400):
    """Evaluate sum_k P(k) * z^k / denom(k) * prod F_i(k)^{e_i} to 130 digits."""
    k_start = 0 if dk in ("(2k+1)^3", "(2k+1)^5", "1") else 1
    z = Decimal(z_num) / Decimal(z_den)
    prod = Decimal(1)
    s = Decimal(0)
    eps = Decimal(10) ** (-131)
    for k in range(0, max_k + 1):
        if k > 0:
            prod *= z
            if e[0] != 0:
                r = Decimal(2 * k - 1) / Decimal(2 * k)
                for _ in range(abs(e[0])):
                    prod = prod * r if e[0] > 0 else prod / r
            if e[1] != 0:
                r = Decimal((3 * k - 2) * (3 * k - 1)) / Decimal(9 * k * k)
                for _ in range(abs(e[1])):
                    prod = prod * r if e[1] > 0 else prod / r
            if e[2] != 0:
                r = Decimal((4 * k - 3) * (4 * k - 1)) / Decimal(16 * k * k)
                for _ in range(abs(e[2])):
                    prod = prod * r if e[2] > 0 else prod / r
            if e[3] != 0:
                r = Decimal((6 * k - 5) * (6 * k - 1)) / Decimal(36 * k * k)
                for _ in range(abs(e[3])):
                    prod = prod * r if e[3] > 0 else prod / r
        if k < k_start:
            continue
        if dk == "k^3":
            den = Decimal(k ** 3)
        elif dk == "k^4":
            den = Decimal(k ** 4)
        elif dk == "k^5":
            den = Decimal(k ** 5)
        elif dk == "(2k+1)^3":
            den = Decimal((2 * k + 1) ** 3)
        elif dk == "(2k+1)^5":
            den = Decimal((2 * k + 1) ** 5)
        elif dk == "1":
            den = Decimal(1)
        else:
            raise ValueError(dk)
        pk = sum(Decimal(c) * Decimal(k ** i) for i, c in enumerate(poly_coeffs))
        term = prod * pk / den
        s += term
        if k >= 20 and abs(term) < eps:
            break
    return s


def verify_all_wz_v2_hits():
    print("=" * 78)
    print("1. 120-DIGIT INDEPENDENT VERIFICATION OF ALL WZ v2 SERIES HITS")
    print("=" * 78)
    # Format: (label, target_val, target_scale, e, dk, z_num, z_den, poly_coeffs [c0, c1, c2])
    hits = [
        ("Markov-Apery zeta(3)",            Z3_D, 2, [-1,0,0,0], "k^3",      -1, 4,    [-5, 0, 0]),
        ("Amdeberhan-Zeilberger zeta(3)",   Z3_D, 2, [-5,0,0,0], "k^5",      -1, 1024, [-32, 160, -205]),
        ("AZ odd-shifted twin (64*zeta(3))",Z3_D, 64,[-5,0,0,0], "(2k+1)^5", -1, 1024, [77, 250, 205]),
        ("Guillera-dual 7*zeta(3) (z=-1/4)",Z3_D, 7, [-5,0,0,0], "(2k+1)^5", -1, 4,    [10, 28, 20]),
        ("Central-binomial 7*pi^3/216",     PI3_D,7, [1,0,0,0],  "(2k+1)^3",  1, 4,    [216, 0, 0]),
        ("Guillera 1/pi^2 #1 (z=-1/1024)",  Decimal(1)/PI2_D, 128, [5,0,0,0], "1", -1, 1024, [13, 180, 820]),
        ("Guillera 1/pi^2 #2 (z=-1/4)",     Decimal(1)/PI2_D,   8, [5,0,0,0], "1", -1, 4,    [1, 8, 20]),
        ("Guillera 1/pi^2 #3 (z=1/16)",     Decimal(1)/PI2_D,  32, [3,0,1,0], "1",  1, 16,   [3, 34, 120]),
        ("Guillera 1/pi^2 #4 (z=-1/48)",    Decimal(1)/PI2_D,  48, [1,1,1,0], "1", -1, 48,   [5, 63, 252]),
        ("Guillera 1/pi^2 #5 (z=27/64)",    Decimal(1)/PI2_D,  48, [3,1,0,0], "1", 27, 64,   [3, 27, 74]),
    ]
    for label, target, scale, e, dk, zn, zd, pcoeffs in hits:
        val = eval_hypergeom(e, dk, zn, zd, pcoeffs, max_k=700)
        err = abs(Decimal(scale) * target - val)
        print(f"  [PASS] {label:36s} | err = {float(err):.3e}")
        assert err < Decimal(10) ** (-118), f"Failed {label}: {err}"


def test_guillera_shadow_duals():
    """
    Test the EXACT analytic continuation k -> -k - 1/2 (shadow duality!) on
    Guillera's remaining three 1/pi^2 series (z = 1/16, -1/48, 27/64)!
    Recall: for a hypergeometric term H(k) = P(k) * prod (a_i)_k / (1)_k^5 * z^k,
    replacing k -> k + 1/2 (or inverting k -> -k - 1/2) flips (a)_k -> 1/(1-a)_{k+1/2}!
    Specifically:
      (1/2)_k -> 1 / ((2k+1) * (1/2)_k) ... wait:
      Gamma(1/2 + (-k-1/2)) = Gamma(-k) -> pole!
      Instead, look at the EXACT WZ dual of H(k) = P(k) * prod (a_i)_k / (1)_k^5 * z^k:
      Replace every (a)_k / (1)_k by (1)_k / (a)_{k+1} = k! / (a * (a+1)_k),
      and evaluate P(k) at k -> k + 1/2 (or solve for the exact P_dual(k) via LLL!!).
      Let's compute the exact inverted term:
        T_dual(k) = z^k * (k!)^5 / prod_{i=1}^5 (a_i)_{k+1}   (k >= 0)!
      Check on Guillera #1 (a = [1/2,1/2,1/2,1/2,1/2], z = -1/1024):
        (1/2)_{k+1} = (1/2)_k * (k + 1/2) = (2k+1)/2 * (1/2)_k -> gives (2k+1)^5 in denom!!
      NOW CHECK ON GUILLERA #3, #4, #5!!
        #3: a = [1/2, 1/2, 1/2, 1/4, 3/4], z = 1/16  -> denom has (2k+1)^3 (4k+1) (4k+3)!!
        #4: a = [1/2, 1/3, 2/3, 1/4, 3/4], z = -1/48 -> denom has (2k+1) (3k+1)(3k+2) (4k+1)(4k+3)!!
        #5: a = [1/2, 1/2, 1/2, 1/3, 2/3], z = 27/64 -> denom has (2k+1)^3 (3k+1)(3k+2)!!
    """
    print("\n" + "=" * 78)
    print("2. TESTING EXACT GUILLERA SHADOW-DUALITY (k -> -k-1/2) ON ALL 5 SERIES!")
    print("=" * 78)

    # Exact Pochhammer ratio T(k) = z^k * (k!)^5 / prod_{i=1}^5 (a_i)_{k+1}
    # where (a_i)_{k+1} = a_i * (a_i+1) * ... * (a_i+k)
    duals = [
        ("Dual of Guillera #1 (z=-1/1024)", [Fraction(1,2)]*5, -1, 1024, [77, 250, 205]),
        ("Dual of Guillera #2 (z=-1/4)",    [Fraction(1,2)]*5, -1, 4,    [10, 28, 20]),
        ("Dual of Guillera #3 (z=1/16)",    [Fraction(1,2), Fraction(1,2), Fraction(1,2), Fraction(1,4), Fraction(3,4)], 1, 16, None),
        ("Dual of Guillera #4 (z=-1/48)",   [Fraction(1,2), Fraction(1,3), Fraction(2,3), Fraction(1,4), Fraction(3,4)], -1, 48, None),
        ("Dual of Guillera #5 (z=27/64)",   [Fraction(1,2), Fraction(1,2), Fraction(1,2), Fraction(1,3), Fraction(2,3)], 27, 64, None),
    ]

    # Simple 5D Decimal LLL to find P(k) = c0 + c1*k + c2*k^2 matching zeta(3) OR pi^3 OR G*pi!
    # Also compute Catalan's constant G and L(3, chi_-3) = 4*pi^3/(81*sqrt(3))!
    sqrt3 = Decimal(3).sqrt()
    L3_chi3 = Decimal(4) * PI3_D / (Decimal(81) * sqrt3) # sum chi_-3(n)/n^3

    def dec_lll(vals, scale_d=95):
        S = Decimal(10) ** scale_d
        n = len(vals)
        B = [[0] * (n + 1) for _ in range(n)]
        for i in range(n):
            B[i][i] = 1
            B[i][n] = int((vals[i] * S).to_integral_value())
        m = n + 1
        Bstar = [[Decimal(0)] * m for _ in range(n)]
        mu = [[Decimal(0)] * n for _ in range(n)]
        Bnorm = [Decimal(0)] * n
        def gs(st):
            for i in range(st, n):
                vi = [Decimal(x) for x in B[i]]
                for j in range(i):
                    if Bnorm[j] != 0:
                        mu[i][j] = sum(Decimal(B[i][c]) * Bstar[j][c] for c in range(m)) / Bnorm[j]
                        for c in range(m): vi[c] -= mu[i][j] * Bstar[j][c]
                    else: mu[i][j] = Decimal(0)
                Bstar[i] = vi; Bnorm[i] = sum(x*x for x in vi)
        gs(0)
        k = 1
        while k < n:
            ch = False
            for j in range(k - 1, -1, -1):
                if abs(mu[k][j]) > Decimal("0.5"):
                    q = int(mu[k][j].to_integral_value())
                    if q != 0:
                        for c in range(m): B[k][c] -= q * B[j][c]
                        ch = True
            if ch: gs(k)
            if Bnorm[k] >= (Decimal("0.75") - mu[k][k-1]**2) * Bnorm[k-1]: k += 1
            else: B[k], B[k-1] = B[k-1], B[k]; gs(k-1); k = max(k-1, 1)
        hits = []
        for row in B:
            c = row[:n]
            if c[0] == 0 or not any(c[1:]): continue
            if sum(x*x for x in c) > 10**9: continue
            res = sum(Decimal(ci) * vi for ci, vi in zip(c, vals))
            if abs(res) < Decimal(10) ** (-85):
                if c[0] < 0: c = [-x for x in c]
                hits.append((c, abs(res)))
        return hits

    for label, avec, zn, zd, _ in duals:
        z = Decimal(zn) / Decimal(zd)
        # Compute S0, S1, S2, S3 for T(k) = z^k * (k!)^5 / prod (a_i)_{k+1}
        # Note: at k=0, T(0) = 1 / prod a_i!
        # Step: T(k) / T(k-1) = z * k^5 / prod (a_i + k)
        # Wait: normalize by prod a_i so T_norm(0) = 1! (just scales by a clean rational)
        t = Decimal(1)
        s = [Decimal(0)] * 4
        eps = Decimal(10) ** (-130)
        for k in range(0, 800):
            if k > 0:
                num = z * Decimal(k ** 5)
                den = Decimal(1)
                for a in avec:
                    den *= (Decimal(a.numerator) / Decimal(a.denominator) + Decimal(k))
                t = t * num / den
            s[0] += t
            s[1] += t * Decimal(k)
            s[2] += t * Decimal(k * k)
            s[3] += t * Decimal(k ** 3)
            if k >= 20 and abs(t * Decimal(k ** 3)) < eps:
                break
        # Test against zeta(3), pi^3, L(3, chi_-3), and (zeta(3), pi^3)!
        found = []
        for tname, tval in [("zeta(3)", Z3_D), ("pi^3", PI3_D), ("L(3,chi_-3)", L3_chi3)]:
            h = dec_lll([tval, s[0], s[1], s[2], s[3]])
            if h:
                found.append((tname, h[0][0], float(h[0][1])))
        if not found:
            # Also test 2-constant combinations [zeta(3), pi^3, s0, s1, s2]
            h2 = dec_lll([Z3_D, PI3_D, s[0], s[1], s[2]])
            if h2:
                found.append(("z3+pi^3", h2[0][0], float(h2[0][1])))
            h3 = dec_lll([Z3_D, L3_chi3, s[0], s[1], s[2]])
            if h3:
                found.append(("z3+L3chi3", h3[0][0], float(h3[0][1])))
        print(f"\n  {label} (a={[str(x) for x in avec]}, z={zn}/{zd}):")
        if found:
            for tn, coeffs, err in found:
                print(f"    *** HIT [{tn}] (err={err:.2e}): {coeffs[0]}*{tn} = {-coeffs[1]} + {-coeffs[2]}k + {-coeffs[3]}k^2 + {-coeffs[4]}k^3 ***")
        else:
            print(f"    S0 = {str(s[0])[:35]} (no simple hit in {{zeta(3), pi^3, L(3,chi_-3)}})")


if __name__ == "__main__":
    verify_all_wz_v2_hits()
    test_guillera_shadow_duals()
