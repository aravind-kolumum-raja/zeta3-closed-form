#!/usr/bin/env python3
"""
check_exp25_joint_exact.py
==========================
Check exact Farey (3 primes/cell) joint mu(A*zeta(3) + B*zeta(4)) on the
two hot candidates surfaced by Exp 25 vs our F_n*(t):
  1. F_n*(t):        R=[(0,13),(2,15),(4,17)], P=(2,23;3,21)  [known: 1st=119.67, full=112.3]
  2. Exp25-Cand-18:  R=[(0,18),(2,20),(4,22)], P=(2,31;4,29)
  3. LK (2610.27799):R=[(0,10),(1,11),(2,12)], P=(1,17;2,16)
  4. Plus nearby L=10..24 z=(0,1,2) and z=(0,2,4) single-presentation candidates!
"""

import math
import cmath
from fractions import Fraction

PI = math.pi


def is_prime(n: int) -> bool:
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True


PRIMES = [p for p in range(43, 2500) if is_prime(p)]


def digamma_p1(z: float) -> float:
    s = 0.0
    for k in range(1, 21): s -= 1.0 / (z + k)
    w = z + 21.0
    return s + math.log(w) - 0.5 / w - 1.0 / (12.0 * w * w) + 1.0 / (120.0 * w ** 4)


def poly_mul(a, b):
    r = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r


def durand_kerner(coeffs):
    c = list(coeffs)
    while len(c) > 1 and abs(c[-1]) < 1e-11: c.pop()
    d = len(c) - 1
    lead = c[d]; c = [x / lead for x in c]
    rad = min(1.0 + max(abs(x) for x in c[:d]), 100.0)
    roots = [rad * cmath.exp(2j * PI * (k + 0.13719) / d) for k in range(d)]
    for _ in range(300):
        mx = 0.0
        for i in range(d):
            z = roots[i]; val = c[d]
            for k in range(d - 1, -1, -1): val = val * z + c[k]
            den = 1.0 + 0j
            for j in range(d):
                if i != j: den *= (z - roots[j])
            if abs(den) < 1e-24: continue
            dz = val / den; roots[i] = z - dz; mx = max(mx, abs(dz))
        if mx < 1e-12: break
    return roots


def xlogx(z: complex) -> complex:
    return 0j if abs(z) < 1e-14 else z * cmath.log(z)


def eval_saddle(rb, s1, p1, s2, p2):
    e1, e2 = p1 - s2, p2 - s1
    pl = [1.0]
    for z, _ in rb: pl = poly_mul(pl, [-z, 1.0])
    pl = poly_mul(pl, [s1, 1.0]); pl = poly_mul(pl, [s2, 1.0])
    pr = [1.0]
    for _, r in rb: pr = poly_mul(pr, [-r, 1.0])
    pr = poly_mul(pr, [p1, 1.0]); pr = poly_mul(pr, [p2, 1.0])
    roots = durand_kerner([a - b for a, b in zip(pl, pr)])
    c_const = e1 * math.log(e1) + e2 * math.log(e2) - sum((r - z) * math.log(r - z) for z, r in rb)
    sigma = -1e100; err_best = -1e100
    for z in roots:
        if abs(z.imag) < 1e-6 and -p2 < z.real < -s2:
            j = -z.real
            sig = c_const + sum((j+r)*math.log(j+r) - ((j+zu)*math.log(j+zu) if j+zu>0 else 0.0) for zu, r in rb)
            sig -= (j-s1)*math.log(j-s1) + (p1-j)*math.log(p1-j) + (j-s2)*math.log(j-s2) + (p2-j)*math.log(p2-j)
            sigma = max(sigma, sig)
        elif z.imag > 1e-6:
            val = complex(c_const, 0.0) + sum(xlogx(z - zu) - xlogx(z - r) for zu, r in rb)
            val -= xlogx(z + p1) - xlogx(z + s1) + xlogx(z + p2) - xlogx(z + s2)
            err_best = max(err_best, val.real - 2.0 * PI * z.imag)
    return sigma, sigma - err_best


def vp_joint_exact(rb, s1, p1, s2, p2, n, p):
    """Exact v_p of (6 d_{M4*n}^4 C_n Q_n, 6 d_{M4*n}^4 C_n P_{3,n}, 6 d_{M4*n}^4 C_n P_{4,n}) mod p^5."""
    zs = sorted(z for z, _ in rb)
    t0_3 = zs[2]
    m4 = p1 + t0_3
    m4n = m4 * n
    max_r = max(r for _, r in rb)
    lim = max(p1 + max_r, m4) * n + 5
    if lim >= p * p: return -1
    K = 5 if p < 600 else 4
    pk = p ** K
    vf = [0]*(lim+1); uf = [1]*(lim+1)
    pH1 = [0]*(lim+1); p2H2 = [0]*(m4n+1); p3H3 = [0]*(m4n+1); p4H4 = [0]*(m4n+1)
    p2_v = p*p; p3_v = p2_v*p; p4_v = p3_v*p
    for i in range(1, lim + 1):
        x, c = i, 0
        while x % p == 0: c += 1; x //= p
        vf[i] = vf[i-1] + c; uf[i] = (uf[i-1]*x) % pk
        inv1 = pow(x, -1, pk)
        pH1[i] = (pH1[i-1] + (inv1 if c > 0 else (p*inv1)%pk)) % pk
        if i <= m4n:
            inv2 = (inv1*inv1)%pk; inv3 = (inv2*inv1)%pk; inv4 = (inv3*inv1)%pk
            p2H2[i] = (p2H2[i-1] + (inv2 if c > 0 else (p2_v*inv2)%pk)) % pk
            p3H3[i] = (p3H3[i-1] + (inv3 if c > 0 else (p3_v*inv3)%pk)) % pk
            p4H4[i] = (p4H4[i-1] + (inv4 if c > 0 else (p4_v*inv4)%pk)) % pk

    e1n = (p1 - s2)*n; e2n = (p2 - s1)*n
    vcn = vf[e1n] + vf[e2n] - sum(vf[(r-z)*n] for z, r in rb)
    ucn_t = (uf[e1n] * uf[e2n]) % pk
    ucn_b = 1
    for z, r in rb: ucn_b = (ucn_b * uf[(r-z)*n]) % pk
    ucn = (ucn_t * pow(ucn_b, -1, pk)) % pk

    inv2 = pow(2, -1, pk); inv3 = pow(3, -1, pk)
    d4_pow = 4 if m4n >= p else 0
    s1n, p1n, s2n, p2n, t0n = s1*n, p1*n, s2*n, p2*n, t0_3*n
    q_ts, p3_ts, p4_ts = [], [], []
    for j in range(s1n, p1n + 1):
        in2 = (s2n <= j <= p2n)
        vn = vcn; ut = ucn; ub = 1
        for zu, ru in rb:
            a, b = j + ru*n, j + zu*n
            vn += vf[a] - vf[b]; ut = (ut*uf[a])%pk; ub = (ub*uf[b])%pk
        a1, b1 = j - s1n, p1n - j
        vn -= vf[a1] + vf[b1]; ub = (ub * uf[a1] % pk * uf[b1]) % pk
        hj = j + t0n
        if in2:
            a2, b2 = j - s2n, p2n - j
            vn -= vf[a2] + vf[b2]; ub = (ub * uf[a2] % pk * uf[b2]) % pk
            uaj = (ut * pow(ub, -1, pk)) % pk
            q_ts.append((vn + d4_pow, uaj))
            # p*(b_j/a_j)
            pld = (sum(pH1[j + zu*n] - pH1[j + ru*n] for zu, ru in rb) + pH1[a1] - pH1[b1] + pH1[a2] - pH1[b2]) % pk
            # P3: a_j H^{(3)} + (1/2) b_j H^{(2)}
            p3_ts.append((vn + d4_pow - 3, (uaj * p3H3[hj]) % pk))
            p3_ts.append((vn + d4_pow - 3, (uaj * pld % pk * p2H2[hj] % pk * inv2) % pk))
            # P4: a_j H^{(4)} + (1/3) b_j H^{(3)}
            p4_ts.append((vn + d4_pow - 4, (uaj * p4H4[hj]) % pk))
            p4_ts.append((vn + d4_pow - 4, (uaj * pld % pk * p3H3[hj] % pk * inv3) % pk))
        else:
            if j < s2n:
                top2, bot2 = p2n - j, s2n - j - 1
                sgn = -1 if ((s2n - j) % 2 == 1) else 1
            else:
                top2, bot2 = j - s2n, j - p2n - 1
                sgn = -1 if ((j - p2n + 1) % 2 == 1) else 1
            vn -= vf[top2] - vf[bot2]; ut = (ut*uf[bot2])%pk; ub = (ub*uf[top2])%pk
            ubj = (ut * pow(ub, -1, pk)) % pk
            if sgn < 0: ubj = (pk - ubj) % pk
            p3_ts.append((vn + d4_pow - 2, (ubj * p2H2[hj] % pk * inv2) % pk))
            p4_ts.append((vn + d4_pow - 3, (ubj * p3H3[hj] % pk * inv3) % pk))

    def eval_vp(ts):
        ts = [(v, u) for v, u in ts if u != 0]
        if not ts: return 99
        vm = min(v for v, _ in ts)
        acc = sum(u * (p ** (v - vm)) for v, u in ts if v - vm < K) % pk
        if acc == 0: return vm + K
        ext = 0
        while acc % p == 0 and ext < K: ext += 1; acc //= p
        return vm + ext

    return min(eval_vp(q_ts), eval_vp(p3_ts), eval_vp(p4_ts))

def cert_joint(name, rb, s1, p1, s2, p2):
    sig, g = eval_saddle(rb, s1, p1, s2, p2)
    zs = sorted(z for z, _ in rb)
    t0_3 = zs[2]; m4 = p1 + t0_3
    max_r = max(r for _, r in rb)
    ma1 = max(p1 + max_r, m4)
    # Farey partition of order max(ma1, m4) + 2
    N = min(max(ma1, m4) + 2, 55)
    fpts = [(0, 1), (1, 1)]
    for d in range(2, N + 1):
        for k in range(1, d):
            if math.gcd(k, d) == 1: fpts.append((k, d))
    fpts.sort(key=lambda x: x[0] / x[1])
    inv_m = 1.0 / m4
    d1 = 0.0
    for i in range(len(fpts) - 1):
        a = fpts[i][0] / fpts[i][1]; b = fpts[i+1][0] / fpts[i+1][1]
        if b <= inv_m + 1e-12: continue
        a_cl = max(a, inv_m)
        if b - a_cl < 1e-12: continue
        min_p = max(int(math.ceil(ma1 * b)) + 6, m4 + 6)
        vals = []
        for p in PRIMES:
            if p < min_p: continue
            mrg = (b - a_cl) * 0.08
            lo = math.ceil((a_cl + mrg) * p); hi = math.floor((b - mrg) * p)
            if lo <= hi and lo >= 1:
                n = (lo + hi) // 2
                if p * p > ma1 * n + 5:
                    v = vp_joint_exact(rb, s1, p1, s2, p2, n, p)
                    if v >= 0: vals.append(v)
                    if len(vals) >= 2: break
        if vals:
            vm = min(vals)
            if vm > 0: d1 += vm * (1.0 / a_cl - 1.0 / b)
    # Periodic u = 1 + x (16 bins for speed)
    p_per = next(p for p in PRIMES if p >= 2 * ma1 + 11)
    d2 = 0.0
    for bi in range(24):
        a, b = bi / 24.0, (bi + 1) / 24.0
        k = max(1, min(p_per - 1, round(0.5 * (a + b) * p_per)))
        v = vp_joint_exact(rb, s1, p1, s2, p2, p_per + k, p_per)
        if v > 0: d2 += v * (digamma_p1(b) - digamma_p1(a))
    dt = d1 + d2
    h4 = sig + 4 * m4 - dt
    surplus = g - h4
    mu4 = g / surplus if surplus > 0 else 999.0
    print(f"  {name:28s}: mu_joint={mu4:8.3f} (surplus={surplus:+7.3f}) | g={g:.3f} sig={sig:.3f} 4M4={4*m4} D4={dt:.4f} ({d1:.3f}+{d2:.3f})")

print("=== EXACT FAREY JOINT mu(A*zeta(3) + B*zeta(4)) CERTIFICATION ===")
cert_joint("F_n*(t) [our paper]", [(0,13),(2,15),(4,17)], 2, 23, 3, 21)
cert_joint("LK (arXiv:2610.27799)", [(0,10),(1,11),(2,12)], 1, 17, 2, 16)
cert_joint("Exp25-Cand-18 (0,2,4)", [(0,18),(2,20),(4,22)], 2, 31, 4, 29)
cert_joint("L=15 z=(0,1,2) P=(1,25;2,23)", [(0,15),(1,16),(2,17)], 1, 25, 2, 23)
cert_joint("L=14 z=(0,1,2) P=(1,23;2,22)", [(0,14),(1,15),(2,16)], 1, 23, 2, 22)
cert_joint("L=12 z=(0,1,2) P=(1,20;2,19)", [(0,12),(1,13),(2,14)], 1, 20, 2, 19)
cert_joint("L=16 z=(0,2,4) P=(2,28;3,25)", [(0,16),(2,18),(4,20)], 2, 28, 3, 25)
