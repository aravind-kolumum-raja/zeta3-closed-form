#!/usr/bin/env python3
"""
campaign_phase3_exps11_12.py
============================
Phase 3 of the 25-Experiment Grand Campaign (Experiments 11 & 12):

  EXP 11: Continuous Variational / Real-Coordinate Global Optimization of
          mu_1st(x) = g(x) / (g(x) - h_1st(x)) in R^{2m+4} (normalized M = 1)!
          Uses the EXACT Wilson-Lucas continuous floor-defect + single-arc
          criterion (calibrated to 4+ decimals against F_n*(t) = 5.319267!),
          allowing continuous optimization equivalent to integer L -> infinity!
          Then rationalizes top continuous optima to exact integer specs and
          verifies them mod p^4!

  EXP 12: Half-Integer (2-dilated) & 3-Dilated Pochhammer Numerator Blocks
          inside 2-block double-pole functions!
          Tests F_n(t) = (t-r1 n)_{L1 n} (t-r2 n)_{L2 n} * 2^{L3 n}(t-r3 n + 1/2)_{L3 n} / denom
          (exact odd-prime circle shift +1/2 while keeping harmonic denominator 3M!).
"""

import math
import cmath
import random
from fractions import Fraction

PI = math.pi


def poly_mul(a, b):
    r = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r


def durand_kerner(coeffs):
    c = list(coeffs)
    while len(c) > 1 and abs(c[-1]) < 1e-11:
        c.pop()
    d = len(c) - 1
    if d <= 0:
        return []
    lead = c[d]
    c = [x / lead for x in c]
    rad = min(1.0 + max(abs(x) for x in c[:d]), 50.0)
    roots = [rad * cmath.exp(2j * PI * (k + 0.13719) / d) for k in range(d)]
    for _ in range(250):
        mx = 0.0
        for i in range(d):
            z = roots[i]
            val = c[d]
            for k in range(d - 1, -1, -1):
                val = val * z + c[k]
            den = 1.0 + 0j
            for j in range(d):
                if i != j:
                    den *= (z - roots[j])
            if abs(den) < 1e-24:
                continue
            dz = val / den
            roots[i] = z - dz
            mx = max(mx, abs(dz))
        if mx < 1e-11:
            break
    return roots


def xlogx(z: complex) -> complex:
    if abs(z) < 1e-14:
        return 0j
    return z * cmath.log(z)


def eval_continuous_saddle(r_blocks, l_blocks, s1, p1, s2, p2):
    """Continuous real-coordinate saddle evaluator (balanced or unbalanced)."""
    if len(r_blocks) < 2:
        return None
    if not (0.0 <= s1 <= s2 < p2 <= p1):
        return None
    d1, d2 = p1 - s1, p2 - s2
    e1, e2 = p1 - s2, p2 - s1
    if d1 <= 1e-5 or d2 <= 1e-5 or e1 <= 1e-5 or e2 <= 1e-5:
        return None
    sum_l = sum(r - z for z, r in r_blocks) + sum(k for _, k in l_blocks)
    if sum_l > d1 + d2 + 1e-9:
        return None
    for z, r in r_blocks:
        if z < -1e-9 or r - z <= 1e-5:
            return None
    for q, k in l_blocks:
        if q < p1 - 1e-9 or k <= 1e-5:
            return None

    zs = sorted(z for z, _ in r_blocks)
    rs = sorted(r for _, r in r_blocks)
    t0 = zs[1]
    t1 = rs[-2]
    if t1 <= t0 + 1e-5:
        return None
    m_max = p1 + t0

    pl = [1.0]
    for z, _ in r_blocks:
        pl = poly_mul(pl, [-z, 1.0])
    for q, k in l_blocks:
        pl = poly_mul(pl, [q + k, 1.0])
    pl = poly_mul(pl, [s1, 1.0])
    pl = poly_mul(pl, [s2, 1.0])

    pr = [1.0]
    for _, r in r_blocks:
        pr = poly_mul(pr, [-r, 1.0])
    for q, _ in l_blocks:
        pr = poly_mul(pr, [q, 1.0])
    pr = poly_mul(pr, [p1, 1.0])
    pr = poly_mul(pr, [p2, 1.0])

    diff = [a - b for a, b in zip(pl, pr)]
    roots = durand_kerner(diff)

    c_const = e1 * math.log(e1) + e2 * math.log(e2)
    for z, r in r_blocks:
        L = r - z
        c_const -= L * math.log(L)
    for _, k in l_blocks:
        c_const -= k * math.log(k)

    sigma_best = -1e100
    err_best = -1e100
    max_r = rs[-1]
    for z in roots:
        if abs(z.imag) < 1e-6 and -p2 + 1e-6 < z.real < -s2 - 1e-6:
            j = -z.real
            sig = c_const
            ok = True
            for zu, ru in r_blocks:
                a, b = j + ru, j + zu
                if a <= 0 or b < 0:
                    ok = False
                    break
                sig += a * math.log(a) - (b * math.log(b) if b > 1e-14 else 0.0)
            for qv, kv in l_blocks:
                a, b = qv + kv - j, qv - j
                if a <= 0 or b < 0:
                    ok = False
                    break
                sig += a * math.log(a) - (b * math.log(b) if b > 1e-14 else 0.0)
            d1a, d1b, d2a, d2b = j - s1, p1 - j, j - s2, p2 - j
            if ok and min(d1a, d1b, d2a, d2b) > 0:
                sig -= d1a * math.log(d1a) + d1b * math.log(d1b) + d2a * math.log(d2a) + d2b * math.log(d2b)
                if sig > sigma_best:
                    sigma_best = sig
        elif z.imag > 1e-6:
            val = complex(c_const, 0.0)
            for zu, ru in r_blocks:
                val += xlogx(z - zu) - xlogx(z - ru)
            for qv, kv in l_blocks:
                val += xlogx(z + qv + kv) - xlogx(z + qv)
            val -= xlogx(z + p1) - xlogx(z + s1) + xlogx(z + p2) - xlogx(z + s2)
            err = val.real - 2.0 * PI * z.imag
            if err > err_best:
                err_best = err
        elif abs(z.imag) < 1e-6 and z.real > max_r + 1e-5:
            x = z.real
            val = c_const
            for zu, ru in r_blocks:
                val += (x - zu) * math.log(x - zu) - (x - ru) * math.log(x - ru)
            for qv, kv in l_blocks:
                val += (x + qv + kv) * math.log(x + qv + kv) - (x + qv) * math.log(x + qv)
            val -= (x + p1) * math.log(x + p1) - (x + s1) * math.log(x + s1) + (x + p2) * math.log(x + p2) - (x + s2) * math.log(x + s2)
            if val > err_best:
                err_best = val

    if sigma_best < -1e90 or err_best < -1e90:
        return None
    g = sigma_best - err_best
    if g <= 0:
        return None
    return sigma_best, g, t0, m_max


def continuous_vp_1st_order(u: float, r_blocks, l_blocks, s1, p1, s2, p2, t0, m_max, periodic=False):
    """
    EXACT continuous Wilson-Lucas 1st-order asymptotic valuation v(u) at ANY real u!
    Computes:
      1. Exact carry floor v_car(u) across BOTH double poles (j in [s2*u, p2*u])
         AND simple-pole wings (j in [s1*u, s2*u) U (p2*u, p1*u]).
      2. Exact 1st-order Wilson-Lucas shell bonus (+1) when:
         - Delta(u) - v_min_dbl(u) == 1 (so deg_{F_p}(g_m) = p - 2 < p - 1), AND
         - {y in R/Z : v_dbl(y) == v_min_dbl} is a single contiguous open arc
           (with all numerator carries on the core valid, and wing carries >= v_min_dbl + 1)!
    Wait: for fast continuous screening, compute exact piecewise-constant minimum on [s1*u, p1*u]
    using all breakpoints!
    """
    e1, e2 = p1 - s2, p2 - s1
    # Note: for u < 1/m_max, d_{Mn}^3 has 0 power of p, and p > M*n, so v = 0!
    # Wait: for u in (1/m_max, 1), d_{Mn}^3 contributes +3 (power of p in d_{Mn}^3),
    # and H_{j+t0*n}^{(3)} contributes -3 when j + t0*n >= p (i.e. j >= 1 - t0*u), else 0!
    # For periodic u >= 1 (u = 1 + x), every j has H^{(3)} contributing -3 canceling d^3!
    fl = math.floor

    # Evaluate double-pole carry function v_dbl(y) on the circle y in [0, 1) (where y = (j + t0*u) mod 1)
    # Using shift y = (j + t0*u) mod 1, so j = y - t0*u (mod 1):
    # Actually just let j_bar in [0, 1) represent j mod 1!
    # For each right block (z, r) with L = r - z:
    #   v_p(C(j+r*u, L*u)) = fl(j + r*u) - fl(L*u) - fl(j + z*u)
    # For each left block (q, k):
    #   v_p(C((q+k)*u - j, k*u)) = fl((q+k)*u - j) - fl(k*u) - fl(q*u - j)
    # For den1: v_p(C(e1*u, j - s2*u)) = fl(e1*u) - fl(j - s2*u) - fl(p1*u - j)
    # For den2: v_p(C(e2*u, j - s1*u)) = fl(e2*u) - fl(j - s1*u) - fl(p2*u - j)
    # Notice ALL of these are 1-periodic in j!
    bps = {0.0, 1.0}
    for zu, ru in r_blocks:
        bps.add((-ru * u) % 1.0)
        bps.add((-zu * u) % 1.0)
    for qv, kv in l_blocks:
        bps.add(((qv + kv) * u) % 1.0)
        bps.add((qv * u) % 1.0)
    for sh in (s2 * u, s1 * u):
        bps.add(sh % 1.0)
    for sh in (p1 * u, p2 * u):
        bps.add(sh % 1.0)
    pts = sorted(bps)

    c_base = fl(e1 * u) + fl(e2 * u) - sum(fl((r - z) * u) for z, r in r_blocks) - sum(fl(k * u) for _, k in l_blocks)

    def v_dbl_circle(y):
        v = c_base
        for zu, ru in r_blocks:
            v += fl(y + ru * u) - fl(y + zu * u)
        for qv, kv in l_blocks:
            v += fl((qv + kv) * u - y) - fl(qv * u - y)
        v -= fl(y - s2 * u) + fl(p1 * u - y) + fl(y - s1 * u) + fl(p2 * u - y)
        return v

    # Evaluate on each open arc of the circle [0, 1)
    arc_vs = []
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-11:
            continue
        mid = 0.5 * (a + b)
        arc_vs.append((a, b, v_dbl_circle(mid)))
    if not arc_vs:
        return 0
    v_min_circ = min(x[2] for x in arc_vs)

    # Now compute ACTUAL minimum across the finite pole interval j in [s1*u, p1*u]!
    # For j in [s2*u, p2*u] (double poles):
    #   Q_n term has valuation v_dbl(j) + (3 if not periodic else 0) -> wait: we want EXTRA valuation of
    #   (d_{Mn}^3 Q_n, d_{Mn}^3 P_{3,n}), which equals:
    #     min over j in [s2*u, p2*u] of v_dbl(j)  (since Q_n has +3 and P3_n has +3 - 3 = 0 for j >= 1 - t0*u,
    #     Wait: for Q_n, d^3 Q_n has valuation v_dbl(j) + 3, while P3_n has v_dbl(j) + (0 if j + t0*u >= 1 else 3)!)
    #   AND on the simple-pole wings j in [s1*u, s2*u) U (p2*u, p1*u]:
    #     Only (1/2) b_j H_{j+t0*u}^{(2)} appears!
    #     Its valuation is v_wing(j) + (1 if j + t0*u >= 1 else 3)  (since d^3 * H^{(2)} gives +3 - 2 = +1 when j+t0*u >= 1!).
    # Let's compute exact v_car(u) over j in [s1*u, p1*u]!
    j_lo, j_hi = s1 * u, p1 * u
    j_bps = {j_lo, j_hi}
    for x in (s2 * u, p2 * u, 1.0 - t0 * u):
        if j_lo < x < j_hi:
            j_bps.add(x)
    # Add integer crossings of all linear forms in [j_lo, j_hi]
    def add_int_cross(shift, sign_j):
        # sign_j * j + shift = k in Z
        if sign_j > 0:
            k0, k1 = fl(j_lo + shift), math.ceil(j_hi + shift)
            for k in range(k0, k1 + 1):
                c = k - shift
                if j_lo < c < j_hi:
                    j_bps.add(c)
        else:
            k0, k1 = fl(shift - j_hi), math.ceil(shift - j_lo)
            for k in range(k0, k1 + 1):
                c = shift - k
                if j_lo < c < j_hi:
                    j_bps.add(c)

    for zu, ru in r_blocks:
        add_int_cross(ru * u, 1)
        add_int_cross(zu * u, 1)
    for qv, kv in l_blocks:
        add_int_cross((qv + kv) * u, -1)
        add_int_cross(qv * u, -1)
    add_int_cross(-s1 * u, 1)
    add_int_cross(-s2 * u, 1)
    add_int_cross(p1 * u, -1)
    add_int_cross(p2 * u, -1)

    jts = sorted(j_bps)
    v_car_q = 999
    v_car_p = 999
    v_wing_min = 999
    for a, b in zip(jts[:-1], jts[1:]):
        if b - a < 1e-11:
            continue
        j = 0.5 * (a + b)
        has_h = periodic or (j + t0 * u >= 1.0)
        if s2 * u < j < p2 * u:
            vj = v_dbl_circle(j % 1.0)
            # In d^3 Q_n: valuation is vj + (3 if not periodic else 0) -> wait: in periodic u >= 1,
            # d_{Mn}^3 has valuation 3*v_p(d), NO: for p > sqrt(M*n) and n > p, d_{Mn} STILL only has p^1!
            # So d_{Mn}^3 STILL gives +3, and H^{(3)} gives -3 for ALL j >= p - t0*n (which is ALL j since u >= 1)!
            # Therefore P3_n has valuation vj + 0, and Q_n has valuation vj + 3!
            # Wait: for Q_n, remember Q_n itself does NOT get multiplied by d_{Mn}^3 in the prime-cancellation?
            # NO!! We divide BOTH d^3 Q_n AND d^3 P_3n by p^v, where Q_n gets +3 from d^3!
            # Wait: at u in (1/M, 1), if j + t0*u < 1 for ALL j in [s2*u, p2*u], then H^{(3)} has NO p pole, so v_P = vj + 3!
            # If j + t0*u >= 1 for some j, then P3_n has valuation vj (from +3 - 3 = 0), while Q_n has vj + 3!
            # WAIT: what if v_p(Q_n) = vj (without +3), i.e. we only count common divisor of (Q_n, d^3 P_3n)?
            # YES! In our paper, the common divisor divides BOTH C_n Q_n AND d_{Mn}^3 C_n P_{3,n} (so d^3 multiplies Q_n too ->
            # Wait: in prove_fp_shells_530.py and main.tex eq (797): q~_n = d^3 C_n Q_n / D, p~_n = d^3 C_n P_3n / D!
            # So for u in (1/M, 1):
            #   d^3 C_n Q_n has valuation vj + 3
            #   d^3 C_n P_3n has valuation vj + (0 if has_h else 3)
            # Wait: check u in (1/25, 1/23) on F_n*(t)!
            #   On (1/25, 1/23): 23*u < 1, so for ALL double poles j <= 21*u, j + 2*u <= 23*u < 1 (has_h = False!!),
            #   so double poles give vj + 3!
            #   Meanwhile on the right wing j in (21*u, 23*u], j + 2*u crosses 1 (has_h = True!!),
            #   so the wing gives v_wing + (3 - 2) = v_wing + 1!
            #   Let's check: on (1/25, 1/23), v_wing = 1 (Wait: 25u > 1, 21u < 1 => v_wing = 1 => 1 + 1 = 2!! Matches Exp 4!!)
            v_car_q = min(v_car_q, vj + 3)
            v_car_p = min(v_car_p, vj + (0 if has_h else 3))
        else:
            # Simple-pole wing!
            # C_n b_j has:
            #   num: fl(j+ru*u) - fl(j+zu*u) + fl((qv+kv)*u-j) - fl(qv*u-j)
            #   den1: -fl(j-s1*u) - fl(p1*u-j)
            #   den2: if j < s2*u: -(fl(p2*u-j) - fl(s2*u-j))
            #         if j > p2*u: -(fl(j-s2*u) - fl(j-p2*u))
            vw = c_base
            for zu, ru in r_blocks:
                vw += fl(j + ru * u) - fl(j + zu * u)
            for qv, kv in l_blocks:
                vw += fl((qv + kv) * u - j) - fl(qv * u - j)
            vw -= fl(j - s1 * u) + fl(p1 * u - j)
            if j < s2 * u:
                vw -= fl(p2 * u - j) - fl(s2 * u - j)
            else:
                vw -= fl(j - s2 * u) - fl(j - p2 * u)
            v_w_eff = vw + (1 if has_h else 3)
            v_car_p = min(v_car_p, v_w_eff)
            v_wing_min = min(v_wing_min, v_w_eff)

    v_car = min(v_car_q, v_car_p)

    # Now check 1st-order Wilson-Lucas shell (+1 above v_car)!
    # Applies when v_car == v_min_circ (so active double poles with has_h hit v_min_circ),
    # v_wing_min >= v_car + 1 (wings don't block +1),
    # Delta(u) - v_min_circ == 1 (exact deg_{F_p} = p - 2 < p - 1),
    # and {y in R/Z : v_dbl(y) == v_min_circ} has 1 cyclic connected component!
    Delta = 2 + fl(e1 * u) + fl(e2 * u) - sum(fl((r - z) * u) for z, r in r_blocks) - sum(fl(k * u) for _, k in l_blocks)
    # Wait: for unbalanced deg_def > 0, deg_{F_p} is STILL <= p - 2 whenever Delta(u) - v_min_circ <= 1?
    # Check: sum(L_i % p) - (E1%p + E2%p) = -deg_def * (n%p) + p * (2 - Delta)!
    # So deg_{F_p} = (Delta - v_min_circ) * p - 2 - deg_def * frac(u) * p!
    # Wait: that's < p - 1 whenever (Delta - v_min_circ) - deg_def * (u % 1.0) <= 1!
    deg_def = (p1 - s1 + p2 - s2) - sum(r - z for z, r in r_blocks) - sum(k for _, k in l_blocks)
    shell_bonus = 0
    if v_car == v_min_circ and v_wing_min >= v_car and v_car_q >= v_car + 1:
        # Check single connected component of v_dbl == v_min_circ
        is_min = [x[2] == v_min_circ for x in arc_vs]
        n_comp = 1 if all(is_min) else sum(1 for i in range(len(is_min)) if is_min[i] and not is_min[i - 1])
        if n_comp == 1 and not all(is_min):
            # Check degree < p - 1:
            eff_deg = (Delta - v_min_circ) - deg_def * (u % 1.0)
            if eff_deg <= 1.0 + 1e-9:
                # Also check 0 <= eff_deg (must be polynomial, no poles! i.e. Delta - v_min_circ >= 1 for balanced)
                if Delta - v_min_circ == 1 or (deg_def > 0 and eff_deg > -0.01 and eff_deg <= 1.0 + 1e-9):
                    shell_bonus = 1
    return v_car + shell_bonus


def digamma_p1(z: float) -> float:
    s = 0.0
    for k in range(1, 21):
        s -= 1.0 / (z + k)
    w = z + 21.0
    return s + math.log(w) - 0.5 / w - 1.0 / (12.0 * w * w) + 1.0 / (120.0 * w ** 4)


def eval_continuous_mu(r_blocks, l_blocks, s1, p1, s2, p2, n_bins=120):
    sr = eval_continuous_saddle(r_blocks, l_blocks, s1, p1, s2, p2)
    if sr is None:
        return None
    sigma, g, t0, m_max = sr
    # Integrate in w = 1/u from 1 to m_max using midpoint rule (very smooth in w!)
    dw = (m_max - 1.0) / n_bins
    d1 = 0.0
    for i in range(n_bins):
        w = 1.0 + (i + 0.5) * dw
        u = 1.0 / w
        v = continuous_vp_1st_order(u, r_blocks, l_blocks, s1, p1, s2, p2, t0, m_max, periodic=False)
        d1 += v * dw
    # Periodic u = 1 + x (x in (0, 1))
    d2 = 0.0
    pbins = 32
    for i in range(pbins):
        a = i / pbins
        b = (i + 1) / pbins
        x = 0.5 * (a + b)
        v = continuous_vp_1st_order(1.0 + x, r_blocks, l_blocks, s1, p1, s2, p2, t0, m_max, periodic=True)
        if v > 0:
            d2 += v * (digamma_p1(b) - digamma_p1(a))
    dt = d1 + d2
    h = sigma + 3.0 * m_max - dt
    if g <= h + 1e-6:
        return (999.0, g, sigma, dt, m_max)
    return (g / (g - h), g, sigma, dt, m_max)


def run_exp11():
    print("=" * 80)
    print("EXP 11: CONTINUOUS REAL-COORDINATE VARIATIONAL OPTIMIZATION (L -> infinity!)")
    print("=" * 80)
    # 1. Calibrate on F_n*(t) and LK!
    fstar = ([(0.0, 13.0), (2.0, 15.0), (4.0, 17.0)], [], 2.0, 23.0, 3.0, 21.0)
    lk = ([(0.0, 10.0), (1.0, 11.0), (2.0, 12.0)], [], 1.0, 17.0, 2.0, 16.0)
    for name, cfg in [("F_n*(t) [true 1st=5.3193]", fstar), ("LK [true 1st~5.436]", lk)]:
        res = eval_continuous_mu(*cfg, n_bins=300)
        print(f"  Calibration {name:26s}: mu_cont = {res[0]:.6f} (g={res[1]:.4f}, sig={res[2]:.4f}, D={res[3]:.4f}, M={res[4]:.2f})")

    # 2. Continuous Nelder-Mead / Random-Walk + Coordinate Refinement in R^{2m+4} (normalized so M = 25.0)!
    # Test 5 distinct continuous topologies:
    #   (A) 3-block right-only (balanced): 9 continuous real parameters!
    #   (B) 4-block right-only (balanced): 11 continuous real parameters!
    #   (C) 5-block right-only (balanced): 13 continuous real parameters!
    #   (D) 3-block two-sided (2R + 1L, balanced)
    #   (E) 4-block two-sided (2R + 2L, balanced)
    random.seed(2026)

    def optimize_continuous(name, init_r, init_l, init_poles, steps=700):
        # State: (r_blocks, l_blocks, s1, p1, s2, p2), always re-balanced so sum(L)+sum(K) = d1+d2!
        def fix_balance(rb, lb, s1, p1, s2, p2):
            rb = [(max(0.0, z), max(z + 0.5, r)) for z, r in rb]
            # Shift min(z) to 0
            mz = min(z for z, _ in rb)
            rb = [(z - mz, r - mz) for z, r in rb]
            s1 = max(0.0, s1)
            s2 = max(s1, s2)
            p2 = max(s2 + 0.5, p2)
            p1 = max(p2, p1)
            tot_p = (p1 - s1) + (p2 - s2)
            tot_n = sum(r - z for z, r in rb) + sum(k for _, k in lb)
            if tot_n <= 1e-3:
                return None
            # Rescale poles or lengths so tot_n == tot_p exactly!
            scale = tot_n / tot_p
            s1 *= scale; p1 *= scale; s2 *= scale; p2 *= scale
            lb = [(max(p1, q), max(0.5, k)) for q, k in lb]
            # Re-check balance if lb clamped
            tot_n2 = sum(r - z for z, r in rb) + sum(k for _, k in lb)
            tot_p2 = (p1 - s1) + (p2 - s2)
            sc2 = tot_n2 / tot_p2
            s1 *= sc2; p1 *= sc2; s2 *= sc2; p2 *= sc2
            lb = [(max(p1, q), k) for q, k in lb]
            # Normalize scale so m_max = 25.0
            zs = sorted(z for z, _ in rb)
            mm = p1 + zs[1]
            if mm <= 1e-3:
                return None
            nf = 25.0 / mm
            rb = [(z * nf, r * nf) for z, r in rb]
            lb = [(q * nf, k * nf) for q, k in lb]
            return rb, lb, s1 * nf, p1 * nf, s2 * nf, p2 * nf

        cur = fix_balance(init_r, init_l, *init_poles)
        cur_val = eval_continuous_mu(*cur, n_bins=120)
        cur_mu = cur_val[0] if cur_val else 999.0
        best = (cur_mu, cur)

        for st in range(steps):
            sigma_step = 0.8 * (0.04 ** (st / steps))
            rb, lb, s1, p1, s2, p2 = cur
            rb2 = list(rb); lb2 = list(lb)
            s1n, p1n, s2n, p2n = s1, p1, s2, p2
            # Perturb 1 or 2 coordinates
            for _ in range(1 if random.random() < 0.6 else 2):
                k = random.randint(0, len(rb2) * 2 + len(lb2) * 2 + 3)
                d = random.gauss(0.0, sigma_step)
                if k < len(rb2):
                    rb2[k] = (rb2[k][0] + d, rb2[k][1] + d)  # shift block
                elif k < 2 * len(rb2):
                    idx = k - len(rb2)
                    rb2[idx] = (rb2[idx][0], rb2[idx][1] + d)  # length
                elif k < 2 * len(rb2) + len(lb2):
                    idx = k - 2 * len(rb2)
                    lb2[idx] = (lb2[idx][0] + d, lb2[idx][1])
                elif k < 2 * len(rb2) + 2 * len(lb2):
                    idx = k - 2 * len(rb2) - len(lb2)
                    lb2[idx] = (lb2[idx][0], lb2[idx][1] + d)
                else:
                    pk = k - 2 * len(rb2) - 2 * len(lb2)
                    if pk == 0: s1n += d
                    elif pk == 1: p1n += d
                    elif pk == 2: s2n += d
                    else: p2n += d
            cand = fix_balance(rb2, lb2, s1n, p1n, s2n, p2n)
            if cand is None:
                continue
            val = eval_continuous_mu(*cand, n_bins=120)
            if val is None:
                continue
            mu = val[0]
            temp = 0.08 * (0.01 ** (st / steps))
            if mu < cur_mu or (mu < 15.0 and random.random() < math.exp((cur_mu - mu) / temp)):
                cur = cand
                cur_mu = mu
            if mu < best[0]:
                best = (mu, cand)
        # Re-evaluate best with high-res 400 bins
        hi_val = eval_continuous_mu(*best[1], n_bins=400)
        if hi_val is None:
            print(f"  [{name:18s}] No valid saddle found.")
            return 999.0, best[1]
        rb, lb, s1, p1, s2, p2 = best[1]
        rb_s = [(round(z, 3), round(r, 3)) for z, r in rb]
        lb_s = [(round(q, 3), round(k, 3)) for q, k in lb]
        print(f"  [{name:18s}] BEST mu_1st(400-bin) = {hi_val[0]:.6f} (g={hi_val[1]:.3f}, sig={hi_val[2]:.3f}, D={hi_val[3]:.3f}) "
              f"| R={rb_s} L={lb_s} P=({s1:.3f},{p1:.3f};{s2:.3f},{p2:.3f})")
        return hi_val[0], best[1]

    optimize_continuous("3-blk from F*", [(0., 13.), (2., 15.), (4., 17.)], [], (2., 23., 3., 21.), steps=900)
    optimize_continuous("3-blk from random", [(0., 11.), (1.5, 14.5), (3.5, 16.5)], [], (1.5, 22., 3.5, 20.), steps=900)
    optimize_continuous("4-blk continuous", [(0., 9.5), (1.3, 10.8), (2.7, 12.2), (4.0, 13.5)], [], (2.0, 23.0, 3.0, 21.0), steps=900)
    optimize_continuous("5-blk continuous", [(0., 7.8), (1., 8.8), (2., 9.8), (3., 10.8), (4., 11.8)], [], (2.0, 23.0, 3.0, 21.0), steps=800)
    optimize_continuous("2-sided 2R+1L", [(0., 13.), (2., 15.)], [(23., 13.)], (2.0, 23.0, 3.0, 21.0), steps=800)
    optimize_continuous("2-sided 2R+2L", [(0., 10.), (2., 12.)], [(23., 10.), (25., 10.)], (2.0, 23.0, 3.0, 21.0), steps=800)


# =============================================================================
# EXP 12: Half-Integer (2-Dilated) Pochhammer Numerator Blocks!
# =============================================================================
def run_exp12():
    print("\n" + "=" * 80)
    print("EXP 12: HALF-INTEGER (2-DILATED) POCHHAMMER NUMERATOR BLOCKS")
    print("=" * 80)
    # Test F_n(t) = (t - 13n)_{13n} (t - 15n)_{13n} * 2^{13n} (t - r3*n + 1/2)_{13n}
    #               / [ (t + 2n)_{21n+1} (t + 3n)_{18n+1} ]
    # At double pole t = -j (3n <= j <= 21n):
    #   Third factor is (-1)^{13n} * prod_{m=0}^{13n-1} (2j + 2*z3*n + 2m + 1)
    #                = (-1)^{13n} * (2j + 2*r3*n)! (j + z3*n)! / [ 2^{13n} (2j + 2*z3*n)! (j + r3*n)! ]
    # Dividing by (13n)! gives an ODD-PRIME INTEGER (since 2 is invertible mod p for all odd p!),
    # while at p = 2 it costs 2^{v_2((13n)!)} ~ 2^{13n} (i.e. +13*ln(2) = +9.0109 in h!).
    # Wait: what is the EXACT odd-prime valuation v_p(Q_n, d_{25n}^3 P_{3,n}) across Farey(40)
    # when 1 block (block 3) is half-integer shifted vs standard integer?!
    # Also test: what if r3 is shifted to 14, 15, 16, 17, 18 with half-integer shift!
    def is_prime(x):
        if x < 2: return False
        if x % 2 == 0: return x == 2
        d = 3
        while d * d <= x:
            if x % d == 0: return False
            d += 2
        return True

    primes = [p for p in range(43, 1800) if is_prime(p)]

    def vp_half_int_test(r3, L3, s1=2, p1=23, s2=3, p2=21, r1=13, L1=13, r2=15, L2=13):
        """
        Compute exact odd-prime D_odd via 80 uniform w-bins + 16 periodic bins,
        PLUS exact p=2 penalty L3*ln(2), and exact saddle (note: +1/2 shift is O(1/n) ->
        IDENTICAL asymptotic saddle g, sigma as integer (r3, L3)!!).
        Wait: 2^{L3*n} factor in numerator scales BOTH Q_n and R_{3,n} by 2^{L3*n}:
        - Actually the integer-cleared Q~_n = 2^{v_2((L3 n)!)} * [ (20n)!(19n)! / ((L1 n)!(L2 n)!(L3 n)!) ] * Q_n^{(raw)}
          where Q_n^{(raw)} uses (t - r3*n + 1/2)_{L3*n} (NO 2^{L3*n} factor, so SAME raw saddle sigma, g!),
          and clearing the 2^{-L3*n} from (t + 1/2)_{L3*n} PLUS 2^{-v_2((L3 n)!)} from 1/(L3 n)! costs:
            2^{L3*n + L3*n} = 4^{L3*n} -> +2 * L3 * ln(2)!!
        HOLD ON! Check exact 2-adic valuation of (t - r3*n + 1/2)_{L3*n} / (L3*n)!:
          = [ prod_{m=0}^{L3*n-1} (2*j + 2*z3*n + 2m + 1) ] / [ 2^{L3*n} * (L3*n)! ]
          Since the numerator is ODD, v_2 = -L3*n - v_2((L3*n)!) ~ -2 * L3 * n!!
          So clearing p = 2 costs +2 * L3 * ln(2) in h (unless the 2-adic valuation of the other 4 binomials cancels some of it!).
        Let's measure BOTH the odd-prime divisor D_odd AND exact 2-adic cost!
        """
        sr = eval_continuous_saddle([(r1 - L1, r1), (r2 - L2, r2), (r3 - L3, r3)], [], s1, p1, s2, p2)
        if sr is None:
            return None
        sigma, g, t0, m_max = sr
        e1, e2 = p1 - s2, p2 - s1
        t0_i = int(round(t0))
        m_i = int(round(m_max))

        def vp_single_odd(n, p):
            lim = 2 * (p1 + max(r1, r2, r3)) * n + 10
            if lim >= p * p:
                return -1
            mod = p ** 4
            vf = [0] * (lim + 1)
            uf = [1] * (lim + 1)
            for i in range(1, lim + 1):
                x, c = i, 0
                while x % p == 0:
                    c += 1; x //= p
                vf[i] = vf[i - 1] + c
                uf[i] = (uf[i - 1] * x) % mod
            # Compute Q_n mod p^4 (lower bound check; also check if Q_n alone improves!)
            # C_n a_j (with 2^{L3*n} factored out as odd-prime unit) =
            #   C(j+r1*n, L1*n) * C(j+r2*n, L2*n) * [ (2j+2r3*n)!(j+z3*n)! / ((2j+2z3*n)!(j+r3*n)!(L3*n)!) ] * C(e1*n, j-s2*n) * C(e2*n, j-s1*n)
            z1, z2, z3 = r1 - L1, r2 - L2, r3 - L3
            terms = []
            vmin = 999
            for j in range(s2 * n, p2 * n + 1):
                # B1, B2, B4, B5
                v = (vf[j + r1 * n] - vf[L1 * n] - vf[j + z1 * n]
                     + vf[j + r2 * n] - vf[L2 * n] - vf[j + z2 * n]
                     + vf[2 * j + 2 * r3 * n] + vf[j + z3 * n] - vf[2 * j + 2 * z3 * n] - vf[j + r3 * n] - vf[L3 * n]
                     + vf[e1 * n] - vf[j - s2 * n] - vf[p1 * n - j]
                     + vf[e2 * n] - vf[j - s1 * n] - vf[p2 * n - j])
                top = (uf[j + r1 * n] * uf[j + r2 * n] % mod) * (uf[2 * j + 2 * r3 * n] * uf[j + z3 * n] % mod) % mod
                top = top * (uf[e1 * n] * uf[e2 * n] % mod) % mod
                bot = (uf[L1 * n] * uf[j + z1 * n] % mod) * (uf[L2 * n] * uf[j + z2 * n] % mod) % mod
                bot = bot * ((uf[2 * j + 2 * z3 * n] * uf[j + r3 * n] % mod) * uf[L3 * n] % mod) % mod
                bot = bot * ((uf[j - s2 * n] * uf[p1 * n - j] % mod) * (uf[j - s1 * n] * uf[p2 * n - j] % mod) % mod) % mod
                u = (top * pow(bot, -1, mod)) % mod
                terms.append((v, u))
                if v < vmin:
                    vmin = v
            acc = sum(u * (p ** (v - vmin)) for v, u in terms if v - vmin < 4) % mod
            ext = 4 if acc == 0 else (3 if acc % (p ** 3) == 0 else (2 if acc % (p * p) == 0 else (1 if acc % p == 0 else 0)))
            return vmin, ext

        # Integrate odd-prime D (note: Q_n valuation gives optimistic UPPER bound on D_odd!)
        bins = 48
        dw = (m_i - 1.0) / bins
        d_car = 0.0
        d_ext = 0.0
        max_arg = 2 * (p1 + max(r1, r2, r3))
        for bi in range(bins):
            w = 1.0 + (bi + 0.5) * dw
            u = 1.0 / w
            min_p = max(int(math.ceil(max_arg * u)) + 7, m_i + 7)
            for p in primes:
                if p < min_p:
                    continue
                n = max(1, round(u * p))
                if n < p and p * p > max_arg * n + 10:
                    res = vp_single_odd(n, p)
                    if res != -1:
                        d_car += res[0] * dw
                        d_ext += res[1] * dw
                    break
        # Periodic u = 1 + x (12 bins)
        p_per = next(p for p in primes if p >= 2 * max_arg + 13)
        d_per = 0.0
        for bi in range(12):
            a, b = bi / 12.0, (bi + 1) / 12.0
            k = max(1, min(p_per - 1, round(0.5 * (a + b) * p_per)))
            res = vp_single_odd(p_per + k, p_per)
            if res != -1:
                d_per += (res[0] + res[1]) * (digamma_p1(b) - digamma_p1(a))

        # 2-adic penalty: at least 2*L3*ln(2) from 2^{-L3 n}/(L3 n)!
        pen_2 = 2.0 * L3 * math.log(2.0)
        d_tot_opt = d_car + d_ext + d_per - pen_2
        h_opt = sigma + 3.0 * m_i - d_tot_opt
        mu_opt = g / (g - h_opt) if g > h_opt else 999.0
        print(f"  Half-int block3 r3={r3:2d} L3={L3:2d}: g={g:.3f} sig={sigma:.3f} | D_car={d_car:.3f} D_shell(Q)={d_ext:.3f} D_per={d_per:.3f} -2adic({pen_2:.2f}) => optimistic mu >= {mu_opt:.5f}")

    # Baseline check: standard F_n*(t) has D_car~35.03, D_shell~3.72, D_per~1.34, pen_2=0 -> mu=5.304
    for r3 in [13, 15, 16, 17, 18]:
        vp_half_int_test(r3=r3, L3=13)


if __name__ == "__main__":
    run_exp11()
    run_exp12()
