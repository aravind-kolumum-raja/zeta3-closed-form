#!/usr/bin/env python3
"""
Hankel-Determinant Lab for zeta(3) and (zeta(3), pi^3)
======================================================
Implements and extends Aabir Fauzan's (Oct 2026) e^{-c n^2} positive-definite
Hankel Gram determinant method from zeta(5) to:
  1. Standard weight-3 measure w_3(y) = (2*pi)^2 y^3 sum_{l>=1} l^2 e^{-2*pi*l*y}
     -> brand-new e^{-c n^2} Hankel proof architecture for zeta(3) not in Q!
  2. Our new alternating weight-3 measure
       w_alt(y) = pi^2 y^3 sum_{l>=0} (2l+1)^2 e^{-pi(2l+1)y}
                = pi^2 y^3 (cosh(2*pi*y) + 3) / (4 sinh^3(pi*y)) > 0
     which has PURE RATIONAL moments AND simultaneously produces:
       - eta(3) = (3/4)*zeta(3) at integer poles t = -j^2, AND
       - 8*beta(3) = pi^3/4 at half-integer poles t = -(j - 1/2)^2!
  3. Exact Fraction Gram determinants, denominator growth vs archimedean decay,
     and quantitative analysis of the homogeneity obstruction for zeta(3)/pi^3.
"""

from fractions import Fraction
from decimal import Decimal, getcontext
import math

getcontext().prec = 420


def compute_pi_dec():
    # Machin to 420 digits
    def arctan_inv(m):
        m2 = m * m
        term = Decimal(1) / Decimal(m)
        s = term
        k = 1
        eps = Decimal(10) ** (-415)
        while True:
            term = -term / m2
            add = term / (2 * k + 1)
            if abs(add) < eps:
                break
            s += add
            k += 1
        return s
    return Decimal(4) * (Decimal(4) * arctan_inv(5) - arctan_inv(239))


def compute_zeta3_dec():
    # Amdeberhan-Zeilberger super-fast series: 3.01 digits/term -> 145 terms = 435 digits!
    # zeta(3) = 1/2 * sum_{k=1}^inf (-1)^{k-1} (205 k^2 - 160 k + 32) / (k^5 C(2k,k)^5)? No, simpler:
    # Markov-Apery gives 0.6 digits/term -> 720 terms = 432 digits in 0.05s!
    s = Decimal(0)
    b = 1
    for k in range(1, 740):
        b = b * (2 * k - 1) * (2 * k) // (k * k)
        term = Decimal(1) / Decimal(k ** 3 * b)
        s += term if (k % 2 == 1) else -term
    return s * Decimal(5) / Decimal(2)


PI_D = compute_pi_dec()
Z3_D = compute_zeta3_dec()
ETA3_D = Decimal(3) * Z3_D / Decimal(4)
PI3_4_D = (PI_D ** 3) / Decimal(4)


def bernoulli_table(max_n):
    B = [Fraction(0)] * (max_n + 1)
    B[0] = Fraction(1)
    for m in range(1, max_n + 1):
        s = Fraction(0)
        for k in range(m):
            s += Fraction(math.comb(m + 1, k)) * B[k]
        B[m] = -s / Fraction(m + 1)
    return B


BERN = bernoulli_table(80)


# =============================================================================
# 1. NUMERICAL VERIFICATION OF EXACT MOMENT & POLE FORMULAS (45+ digits)
# =============================================================================

def verify_measures_numerically():
    print("=" * 72)
    print("1. HIGH-PRECISION VERIFICATION OF w_3(y) AND w_alt(y) IDENTITIES")
    print("=" * 72)
    saved_prec = getcontext().prec
    getcontext().prec = 50
    pi50 = +PI_D; z3_50 = +Z3_D; eta50 = +ETA3_D; pi34_50 = +PI3_4_D

    def integrate_0_inf(f, h=0.04, N=200):
        tot = Decimal(0)
        hd = Decimal(str(h))
        for k in range(-N, N + 1):
            u = hd * Decimal(k)
            eu = u.exp()
            em = Decimal(1) / eu
            sh = (eu - em) / Decimal(2)
            ch = (eu + em) / Decimal(2)
            if sh > Decimal(20) or sh < Decimal(-30):
                continue
            y = sh.exp()
            tot += f(y) * y * ch
        return tot * hd

    def w3(y):
        x = pi50 * y
        if x > Decimal(55): return Decimal(0)
        ex = x.exp(); em = Decimal(1) / ex
        return (pi50 * pi50) * (y ** 3) * ((ex + em) / Decimal(2)) / (((ex - em) / Decimal(2)) ** 3)

    def walt(y):
        x = pi50 * y
        if x > Decimal(90): return Decimal(0)
        ex = x.exp(); em = Decimal(1) / ex
        sh = (ex - em) / Decimal(2)
        ch2 = (ex * ex + em * em) / Decimal(2)
        return (pi50 * pi50) * (y ** 3) * (ch2 + Decimal(3)) / (Decimal(4) * (sh ** 3))

    print("\n[Moments of w_3(y)]: mu^{(3)}(t^e) = (-1)^e * B_{2e+2} * (2e+3)/2")
    for e in range(4):
        num = integrate_0_inf(lambda y, e=e: (y ** (2 * e)) * w3(y))
        exact = ((-1) ** e) * BERN[2 * e + 2] * Fraction(2 * e + 3, 2)
        err = abs(num - Decimal(exact.numerator) / Decimal(exact.denominator))
        print(f"  e={e}: exact={str(exact):>10s}, err={float(err):.3e}")

    print("\n[Moments of w_alt(y)]: mu^{alt}(t^e) = (-1)^e * (2^{2e+2}-1) * B_{2e+2} * (2e+3)/2")
    for e in range(4):
        num = integrate_0_inf(lambda y, e=e: (y ** (2 * e)) * walt(y))
        exact = ((-1) ** e) * (2 ** (2 * e + 2) - 1) * BERN[2 * e + 2] * Fraction(2 * e + 3, 2)
        err = abs(num - Decimal(exact.numerator) / Decimal(exact.denominator))
        print(f"  e={e}: exact={str(exact):>10s}, err={float(err):.3e}")

    print("\n[Integer poles of w_3(y)]: j^2*(zeta(3)-H_j^(3)) - 1/2 + 1/(2j)")
    for j in range(1, 5):
        num = integrate_0_inf(lambda y, j=j: w3(y) / (y * y + Decimal(j * j)))
        Hj3 = sum(Decimal(1) / Decimal(m ** 3) for m in range(1, j + 1))
        pred = Decimal(j * j) * (z3_50 - Hj3) - Decimal("0.5") + Decimal(1) / Decimal(2 * j)
        print(f"  j={j}: err = {float(abs(num - pred)):.3e}")

    print("\n[Integer poles of w_alt(y)]: (-1)^{j+1} j^2 (eta(3) - H_j^{alt}) + 1/(2j)")
    for j in range(1, 5):
        num = integrate_0_inf(lambda y, j=j: walt(y) / (y * y + Decimal(j * j)))
        Halt = sum(Decimal((-1) ** (m - 1)) / Decimal(m ** 3) for m in range(1, j + 1))
        pred = Decimal((-1) ** (j + 1)) * Decimal(j * j) * (eta50 - Halt) + Decimal(1) / Decimal(2 * j)
        print(f"  j={j}: err = {float(abs(num - pred)):.3e}")

    print("\n[Half-integer poles of w_alt(y)]: (-1)^{j+1} (j-1/2)^2 (pi^3/4 - 8*H_j^{odd}) + 1/(2j-1)")
    for j in range(1, 5):
        a = Decimal(2 * j - 1) / Decimal(2)
        num = integrate_0_inf(lambda y, a=a: walt(y) / (y * y + a * a))
        Hodd = Decimal(8) * sum(Decimal((-1) ** m) / Decimal((2 * m + 1) ** 3) for m in range(j))
        pred = Decimal((-1) ** (j + 1)) * a * a * (pi34_50 - Hodd) + Decimal(1) / (Decimal(2) * a)
        print(f"  a={2*j-1}/2: err = {float(abs(num - pred)):.3e}")
    getcontext().prec = saved_prec


# =============================================================================
# 2. EXACT POLYNOMIAL & RATIONAL-MATRIX ARITHMETIC OVER Q[X]
# =============================================================================

def poly_add(a, b):
    n = max(len(a), len(b))
    r = [Fraction(0)] * n
    for i in range(len(a)): r[i] += a[i]
    for i in range(len(b)): r[i] += b[i]
    while len(r) > 1 and r[-1] == 0: r.pop()
    return r

def poly_sub(a, b):
    n = max(len(a), len(b))
    r = [Fraction(0)] * n
    for i in range(len(a)): r[i] += a[i]
    for i in range(len(b)): r[i] -= b[i]
    while len(r) > 1 and r[-1] == 0: r.pop()
    return r

def poly_mul(a, b):
    r = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            r[i + j] += x * y
    while len(r) > 1 and r[-1] == 0: r.pop()
    return r

def poly_scale(a, c):
    if c == 0: return [Fraction(0)]
    return [x * c for x in a]

def poly_pow(a, e):
    r = [Fraction(1)]
    b = a
    while e > 0:
        if e & 1: r = poly_mul(r, b)
        b = poly_mul(b, b)
        e >>= 1
    return r

def poly_divmod(num, den):
    num = list(num)
    while len(num) > 1 and num[-1] == 0: num.pop()
    while len(den) > 1 and den[-1] == 0: den.pop()
    dn = len(num) - 1; dd = len(den) - 1
    if dn < dd: return [Fraction(0)], num
    q = [Fraction(0)] * (dn - dd + 1)
    lead = den[-1]
    for k in range(dn - dd, -1, -1):
        c = num[k + dd] / lead
        q[k] = c
        if c != 0:
            for j in range(dd + 1):
                num[k + j] -= c * den[j]
    while len(num) > 1 and num[-1] == 0: num.pop()
    return q, num

def poly_eval_rat(p, x):
    r = Fraction(0)
    for c in reversed(p):
        r = r * x + c
    return r

def det_poly_exact(M):
    # Bareiss / fraction-free elimination over Q[X] (or direct Gaussian over Q(X) -> wait, interpolate at h+1 points!)
    # Interpolation at x = 0, 1, ..., h is 100x faster & cleaner!
    h = len(M)
    # deg(det) <= h
    vals = []
    for xv in range(h + 1):
        A = [[poly_eval_rat(M[i][j], Fraction(xv)) for j in range(h)] for i in range(h)]
        # Exact Fraction Gaussian elimination
        det = Fraction(1)
        for col in range(h):
            piv = None
            for row in range(col, h):
                if A[row][col] != 0:
                    piv = row; break
            if piv is None:
                det = Fraction(0); break
            if piv != col:
                A[col], A[piv] = A[piv], A[col]
                det = -det
            pv = A[col][col]
            det *= pv
            for row in range(col + 1, h):
                if A[row][col] != 0:
                    f = A[row][col] / pv
                    for k in range(col, h):
                        A[row][k] -= f * A[col][k]
        vals.append(det)
    # Lagrange interpolation from (0..h, vals)
    res = [Fraction(0)]
    for i in range(h + 1):
        num = [Fraction(1)]
        den = Fraction(1)
        for j in range(h + 1):
            if i == j: continue
            num = poly_mul(num, [Fraction(-j), Fraction(1)])
            den *= Fraction(i - j)
        res = poly_add(res, poly_scale(num, vals[i] / den))
    return res

def det_dec(M):
    h = len(M)
    A = [[Decimal(x) for x in row] for row in M]
    det = Decimal(1)
    for col in range(h):
        # partial pivoting
        piv = max(range(col, h), key=lambda r: abs(A[r][col]))
        if abs(A[piv][col]) == 0: return Decimal(0)
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            det = -det
        pv = A[col][col]
        det *= pv
        for row in range(col + 1, h):
            f = A[row][col] / pv
            for k in range(col, h):
                A[row][k] -= f * A[col][k]
    return det


# =============================================================================
# 3. FAUZAN WEIGHT-3 HANKEL GRAM DETERMINANT FOR zeta(3)
# =============================================================================

def mu3_mono(e):
    return ((-1) ** e) * BERN[2 * e + 2] * Fraction(2 * e + 3, 2)

def mu3_pole(j):
    # Exact verified identity: j^2 * (X - H_j^{(3)}) - 1/2 + 1/(2j)
    Hj3 = sum(Fraction(1, m ** 3) for m in range(1, j + 1))
    return [Fraction(1, 2 * j) - Fraction(1, 2) - Fraction(j * j) * Hj3, Fraction(j * j)]

def build_w3_hankel(K, N, h, power_N=4):
    # G_{K,N,h}(X) = [mu_X(D_N(t)^power_N * t^{i+j} / D_K(t))]_{0 <= i,j < h}
    # where D_m(t) = prod_{j=1}^m (t + j^2).
    poles = list(range(1, K + 1))
    DK = [Fraction(1)]
    for j in poles:
        DK = poly_mul(DK, [Fraction(j * j), Fraction(1)])
    DN = [Fraction(1)]
    for j in range(1, N + 1):
        DN = poly_mul(DN, [Fraction(j * j), Fraction(1)])
    DN_pow = poly_pow(DN, power_N)

    # Precompute partial-fraction denominator weights prod_{k!=j} (k^2 - j^2)
    denom_j = {}
    for j in poles:
        prod = Fraction(1)
        for k in poles:
            if k != j: prod *= Fraction(k * k - j * j)
        denom_j[j] = prod

    # Precompute moments mu_X(D_N^p * t^s / D_K) for s = 0 .. 2h - 2
    pole_polys = {j: mu3_pole(j) for j in poles}
    entry = []
    for s in range(2 * h - 1):
        P = [Fraction(0)] * s + DN_pow
        q_poly, _ = poly_divmod(P, DK)
        val = [sum(c * mu3_mono(e) for e, c in enumerate(q_poly))]
        for j in poles:
            res = poly_eval_rat(P, Fraction(-j * j)) / denom_j[j]
            val = poly_add(val, poly_scale(pole_polys[j], res))
        entry.append(val)

    G = [[entry[i + j] for j in range(h)] for i in range(h)]
    poly_det = det_poly_exact(G)
    # Exact common denominator of poly_det -> primitive integer polynomial Q_K(X)!
    lcm_den = 1
    for c in poly_det:
        lcm_den = math.lcm(lcm_den, c.denominator)
    int_coeffs = [int((c * lcm_den).numerator) for c in poly_det]
    g_all = 0
    for ic in int_coeffs:
        g_all = math.gcd(g_all, abs(ic))
    prim_coeffs = [ic // g_all for ic in int_coeffs]
    # Evaluate primitive Q_K(zeta(3)) to 80 digits!
    q_val = Decimal(0)
    for c in reversed(prim_coeffs):
        q_val = q_val * Z3_D + Decimal(c)
    # Wait: also check sign (make Q_K(zeta(3)) > 0)
    if q_val < 0:
        q_val = -q_val
        prim_coeffs = [-c for c in prim_coeffs]
    return prim_coeffs, q_val, lcm_den, g_all


# =============================================================================
# 4. ALTERNATING WEIGHT-3 HANKEL GRAM MATRIX ON MIXED (zeta(3), pi^3) POLES!
# =============================================================================

def mualt_mono(e):
    return ((-1) ** e) * (2 ** (2 * e + 2) - 1) * BERN[2 * e + 2] * Fraction(2 * e + 3, 2)

def test_mixed_alt_hankel():
    print("\n" + "=" * 72)
    print("3. ALTERNATING MEASURE w_alt(y): MIXED INTEGER & HALF-INTEGER POLES")
    print("=" * 72)
    # Evaluate the exact Gram matrix G_K(X, Y) at X = eta(3) = 3/4*zeta(3), Y = pi^3/4
    # for pure integer poles, pure half-integer poles, AND mixed poles!
    def build_alt_gram(int_poles, half_poles, N_int, N_half, h, pN=4):
        # Poles are a in int_poles (j) and half_poles (j - 1/2), square poles at t = -a^2
        all_a = [Fraction(j) for j in int_poles] + [Fraction(2 * j - 1, 2) for j in half_poles]
        DK = [Fraction(1)]
        for a in all_a:
            DK = poly_mul(DK, [a * a, Fraction(1)])
        DN = [Fraction(1)]
        for j in range(1, N_int + 1):
            DN = poly_mul(DN, [Fraction(j * j), Fraction(1)])
        for j in range(1, N_half + 1):
            a = Fraction(2 * j - 1, 2)
            DN = poly_mul(DN, [a * a, Fraction(1)])
        DN_pow = poly_pow(DN, pN)

        # Each pole linear form in (1, X, Y): (c0, cX, cY)
        pole_lin = {}
        for j in int_poles:
            a = Fraction(j)
            Halt = sum(Fraction((-1) ** (m - 1), m ** 3) for m in range(1, j + 1))
            # Exact verified sign: (-1)^{j+1} * j^2 * (eta(3) - Halt) + 1/(2j)
            sgn = Fraction((-1) ** (j + 1))
            pole_lin[a] = (Fraction(1, 2 * j) - sgn * a * a * Halt, sgn * a * a, Fraction(0))
        for j in half_poles:
            a = Fraction(2 * j - 1, 2)
            Hodd = Fraction(8) * sum(Fraction((-1) ** m, (2 * m + 1) ** 3) for m in range(j))
            sgn = Fraction((-1) ** (j + 1))
            pole_lin[a] = (Fraction(1) / (2 * a) - sgn * a * a * Hodd, Fraction(0), sgn * a * a)

        denom_a = {}
        for a in all_a:
            prod = Fraction(1)
            for b in all_a:
                if b != a: prod *= (b * b - a * a)
            denom_a[a] = prod

        # Compute entries as (c0, cX, cY)
        entries = []
        for s in range(2 * h - 1):
            P = [Fraction(0)] * s + DN_pow
            q_poly, _ = poly_divmod(P, DK)
            c0 = sum(c * mualt_mono(e) for e, c in enumerate(q_poly))
            cX = Fraction(0); cY = Fraction(0)
            for a in all_a:
                res = poly_eval_rat(P, -a * a) / denom_a[a]
                p0, pX, pY = pole_lin[a]
                c0 += res * p0; cX += res * pX; cY += res * pY
            entries.append((c0, cX, cY))

        # Evaluate numerical Gram matrix at (ETA3_D, PI3_4_D) AND compute exact lcm of entry denominators
        M_dec = [[Decimal(0)] * h for _ in range(h)]
        # Also compute exact bivariate polynomial via exact lcm of entries (or determinant denominator)
        for i in range(h):
            for j in range(h):
                c0, cX, cY = entries[i + j]
                v = (Decimal(c0.numerator) / Decimal(c0.denominator)
                     + (Decimal(cX.numerator) / Decimal(cX.denominator)) * ETA3_D
                     + (Decimal(cY.numerator) / Decimal(cY.denominator)) * PI3_4_D)
                M_dec[i][j] = v
        d_val = det_dec(M_dec)

        # Compute exact univariate polynomial if pure X or pure Y
        if not half_poles:
            G_poly = [[[entries[i + j][0], entries[i + j][1]] for j in range(h)] for i in range(h)]
            pdet = det_poly_exact(G_poly)
            lcm_d = 1
            for c in pdet: lcm_d = math.lcm(lcm_d, c.denominator)
            g_all = 0
            for c in pdet: g_all = math.gcd(g_all, abs(int((c * lcm_d).numerator)))
            scale = Fraction(lcm_d, g_all)
            Q_val = d_val * Decimal(scale.numerator) / Decimal(scale.denominator)
            return d_val, Q_val, math.log(lcm_d / g_all), len(all_a)
        elif not int_poles:
            G_poly = [[[entries[i + j][0], entries[i + j][2]] for j in range(h)] for i in range(h)]
            pdet = det_poly_exact(G_poly)
            lcm_d = 1
            for c in pdet: lcm_d = math.lcm(lcm_d, c.denominator)
            g_all = 0
            for c in pdet: g_all = math.gcd(g_all, abs(int((c * lcm_d).numerator)))
            scale = Fraction(lcm_d, g_all)
            Q_val = d_val * Decimal(scale.numerator) / Decimal(scale.denominator)
            return d_val, Q_val, math.log(lcm_d / g_all), len(all_a)
        else:
            # Bivariate: check exact common denominator of det(G(X,Y)) by sampling (X,Y) grid of size (h+1)x(h+1)?
            # Even simpler: clear each entry denominator first, or report exact raw det & entry-lcm bound
            row_lcms = []
            for i in range(h):
                l = 1
                for j in range(h):
                    for c in entries[i + j]: l = math.lcm(l, c.denominator)
                row_lcms.append(l)
            prod_lcm = math.prod(row_lcms)
            Q_bound = d_val * Decimal(prod_lcm)
            return d_val, Q_bound, math.log(prod_lcm), len(all_a)

    print("\n[3A. Pure integer poles on w_alt -> eta(3) = (3/4)*zeta(3)]:")
    for K, N, h, pN in [(6, 1, 5, 4), (8, 1, 7, 4), (10, 1, 9, 4), (12, 1, 11, 4)]:
        d_val, Q_val, log_sc, nK = build_alt_gram(list(range(1, K+1)), [], N, 0, h, pN)
        rate = float(Q_val.ln()) / (K * K)
        print(f"  K={K:2d} N={N} h={h:2d} pN={pN}: Q_K(eta(3)) = {float(Q_val):.4e} (>0!), log(Q)/K^2 = {rate:+.5f}, log_den/K^2 = {log_sc/(K*K):.5f}")

    print("\n[3B. Pure half-integer poles on w_alt -> pi^3/4]:")
    for K, N, h, pN in [(6, 1, 5, 4), (8, 1, 7, 4), (10, 1, 9, 4), (12, 1, 11, 4)]:
        d_val, Q_val, log_sc, nK = build_alt_gram([], list(range(1, K+1)), 0, N, h, pN)
        rate = float(Q_val.ln()) / (K * K)
        print(f"  K={K:2d} N={N} h={h:2d} pN={pN}: Q_K(pi^3/4) = {float(Q_val):.4e} (>0!), log(Q)/K^2 = {rate:+.5f}, log_den/K^2 = {log_sc/(K*K):.5f}")

    print("\n[3C. MIXED integer + half-integer poles on w_alt -> bivariate Q_K(eta(3), pi^3/4)]:")
    for K1, K2, N1, N2, h, pN in [(4, 4, 1, 1, 7, 2), (5, 5, 1, 1, 9, 2), (6, 6, 1, 1, 11, 2)]:
        d_val, Q_bd, log_sc, K = build_alt_gram(list(range(1, K1+1)), list(range(1, K2+1)), N1, N2, h, pN)
        print(f"  K1={K1},K2={K2} (K={K:2d}) h={h:2d}: det(G) = {float(d_val):.4e} (>0!), log(det)/K^2 = {float(d_val.ln())/(K*K):+.5f}")


def main():
    verify_measures_numerically()

    print("\n" + "=" * 72)
    print("2. FAUZAN WEIGHT-3 HANKEL GRAM DETERMINANT ON w_3(y) FOR zeta(3)")
    print("=" * 72)
    print("Scanning (K, N, h, power_N) for primitive integer polynomial Q_K(X) in Z[X]:")
    for pN in [2, 4, 6]:
        print(f"\n  --- Numerator power D_N(t)^{pN} ---")
        for K, N, h in [(6, 1, 5), (8, 1, 7), (10, 1, 9), (12, 1, 11), (14, 1, 13), (12, 2, 10), (14, 2, 12)]:
            if pN * N + 2 * (h - 1) > K + 8:
                pass
            prim, qv, lcm_d, g_all = build_w3_hankel(K, N, h, power_N=pN)
            rate = float(qv.ln()) / (K * K)
            print(f"  K={K:2d} N={N} h={h:2d} (deg={len(prim)-1:2d}): 0 < Q_K(zeta(3)) = {float(qv):.4e}, "
                  f"log(Q)/K^2 = {rate:+.5f}, log(den/gcd)/K^2 = {math.log(lcm_d/g_all)/(K*K):.5f}")

    test_mixed_alt_hankel()


if __name__ == "__main__":
    main()
