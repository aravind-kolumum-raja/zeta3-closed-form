#!/usr/bin/env python3
"""
certify_rv_bound.py — Independent 80+ digit exact-rational & rigorous-error-bound
certifier for Theorem 5.1 & Lemma 5.1 of Rhin & Viola, "The group structure for
zeta(3)", Acta Arith. 97 (2001), 269-293.

Uses ONLY Python standard library (fractions, decimal, math) — zero floating-point
cut-point ambiguity (all breakpoints compared as exact integers a*d' vs a'*d),
110-digit Decimal arithmetic for stationary points and digamma psi(x), and the
exact alternating-series remainder bound of the Bernoulli asymptotic expansion
for psi(x) (strictly alternating for x > 0) to produce a rigorous enclosure
[mu_lo, mu_hi].

Usage:
    python3 scripts/certify_rv_bound.py
    python3 scripts/certify_rv_bound.py 1167 1240 1386 1094 875 802 656 948
"""

import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction

# Set 110 decimal digits of working precision
getcontext().prec = 110

# Six families of level-4 cosets of Phi = W(D_5) mod Theta (RV01, Sec. 4)
FN = (
    (0, 2, 3, 6),
    (0, 1, 3, 5),
    (0, 3, 4, 7),
    (1, 2, 5, 6),
    (0, 1, 4, 5),
    (1, 3, 5, 7),
)
FD = (
    (8, 12, 13, 14),
    (11, 13, 14, 15),
    (9, 10, 13, 14),
    (8, 9, 12, 13),
    (9, 10, 13, 14),
    (9, 11, 13, 15),
)
FR = (8, 8, 4, 4, 4, 2)


def theta(i: int, t: int) -> int:
    return (i + t) % 8 if i < 8 else 8 + ((i - 8 + t) % 8)


def bernoulli_even_exact(max_k: int):
    """Exact even Bernoulli numbers B_2, B_4, ..., B_{2*max_k} as Fractions."""
    # Akiyama-Tanigawa algorithm for B_n
    n_max = 2 * max_k
    A = [Fraction(0)] * (n_max + 1)
    B = []
    for m in range(n_max + 1):
        A[m] = Fraction(1, m + 1)
        for j in range(m, 0, -1):
            A[j - 1] = j * (A[j - 1] - A[j])
        B.append(A[0])
    # Return list of (2k, Decimal(B_{2k} / (2k))) for k = 1..max_k
    out = []
    for k in range(1, max_k + 1):
        bk = B[2 * k] / (2 * k)
        out.append((2 * k, Decimal(bk.numerator) / Decimal(bk.denominator)))
    return out


# Precompute B_{2k}/(2k) up to 2k = 60 (30 terms).
# For y >= 45, the 30th term |B_60|/(60 * 45^60) < 3e-70; wait: for y >= 60:
# |B_60| ~ 1.2e30, 60^60 ~ 4.8e106 => term ~ 4e-79!
# Let's use shift M_SHIFT = 64 and K_TERMS = 32:
# |B_64| ~ 1.07e34, 64^64 ~ 3.9e115 => error < 1e-82 per psi evaluation!
M_SHIFT = 64
K_TERMS = 32
BERN_COEFFS = bernoulli_even_exact(K_TERMS + 1)
BERN_USE = BERN_COEFFS[:K_TERMS]
BERN_NEXT_ABS = abs(BERN_COEFFS[K_TERMS][1])  # |B_{2K+2}/(2K+2)| for exact remainder bound


def psi_dec(x: Decimal) -> tuple[Decimal, Decimal]:
    """
    Compute psi(x) = Gamma'(x)/Gamma(x) for x in (0, 1] with a rigorous upper bound
    on absolute error (truncation + rounding).

    By the Euler-Maclaurin / Stieltjes alternating remainder theorem for digamma,
    for y > 0, the remainder after K terms has strictly smaller magnitude than
    the (K+1)-th term |B_{2K+2}| / ((2K+2) * y^{2K+2}).
    """
    acc = Decimal(0)
    for t in range(M_SHIFT):
        acc -= Decimal(1) / (x + t)
    y = x + M_SHIFT
    inv_y = Decimal(1) / y
    inv_y2 = inv_y * inv_y
    out = y.ln() - Decimal("0.5") * inv_y
    p = inv_y2
    for _, coeff in BERN_USE:
        out -= coeff * p
        p *= inv_y2
    # Rigorous truncation bound + generous 1e-95 rounding margin at prec=110
    trunc_err = BERN_NEXT_ABS * p
    err = trunc_err + Decimal("1e-95")
    return out + acc, err


def certify(P: tuple[int, ...], title: str = ""):
    t0 = time.perf_counter()
    assert len(P) == 8, "Need 8 integers (h,j,k,l,m,q,r,s)"
    h, j, k, l, m, q, r, s = P
    print(f"\n================ {title or 'CERTIFICATION'} ================")
    print(f"P = (h,j,k,l,m,q,r,s) = {P}")

    # 1. Exact linear relations (2.2) and (2.3)
    assert h + m == k + r, f"Violated (2.2): h+m={h+m} != k+r={k+r}"
    assert j + q == l + s, f"Violated (2.3): j+q={j+q} != l+s={l+s}"

    # 2. Auxiliary integers (4.7) and BOTH dual expressions
    V = [
        h, j, k, l, m, q, r, s,
        h + l - j, j + m - k, k + q - l, l + r - m,
        m + s - q, q + h - r, r + j - s, s + k - h,
    ]
    dual = [h + q - s, j + r - h, k + s - j, l + h - k, m + j - l, q + k - m, r + l - q, s + m - r]
    assert V[8:] == dual, "Violated dual expressions (4.7)"
    assert all(v > 0 for v in V), f"Non-positive integer in V={V}"
    print(f"V' = (h',j',k',l',m',q',r',s') = {tuple(V[8:])}")

    top = sorted(V, reverse=True)
    M, N, Q = top[0], top[1], top[2]
    print(f"M, N, Q = {M}, {N}, {Q} (sum = {M + N + Q})")

    # 3. Stationary points (x0,y0,z0) and (x1,y1,z1) to 110 digits + gradient verification
    A = Decimal(V[8] * V[11] * V[14])
    B = Decimal((h - j) * (h - k) * (q - r) + (h * l + l * r + r * h) * V[15] - h * l * r)
    C = Decimal(-h * r * V[15])
    disc = (B * B - Decimal(4) * A * C).sqrt()
    dh, dj, dk, dl, dq, ds = map(Decimal, (h, j, k, l, q, s))
    dv8, dv11, dv13, dv15 = map(Decimal, (V[8], V[11], V[13], V[15]))

    c0 = Decimal(0)
    c1 = Decimal(0)
    for e in (0, 1):
        x = (-B + disc) / (Decimal(2) * A) if e == 0 else (-B - disc) / (Decimal(2) * A)
        y = (dv11 * x + (dk - dh)) / ((dv11 - ds) * x + dv15)
        z = (
            dv8 * dv11 * x * x
            - (dh * dv8 + (dh - ds) * dv11 - dk * (dh + dq)) * x
            - dh * dv15
        ) / ((x - Decimal(1)) * ((dh - ds) * dv11 * x + dh * dv15))

        if e == 0:
            assert Decimal(0) < x < Decimal(1), "x0 outside (0,1)"
            assert Decimal(0) < y < Decimal(1), "y0 outside (0,1)"
            assert Decimal(0) < z < Decimal(1), "z0 outside (0,1)"
        else:
            assert x < Decimal(0) and y < Decimal(0), "x1,y1 not < 0"
            assert x * y > Decimal(1), "x1*y1 not > 1"
            assert z < Decimal(1) / (Decimal(1) - x * y), "z1 not < 1/(1-x1*y1)"

        # Verify grad log|f| = 0 directly
        D = Decimal(1) - (Decimal(1) - x * y) * z
        gx = dh / x - dl / (Decimal(1) - x) - dv13 * (y * z) / D
        gy = dk / y - ds / (Decimal(1) - y) - dv13 * (x * z) / D
        gz = dj / z - dq / (Decimal(1) - z) + dv13 * (Decimal(1) - x * y) / D
        grad_max = max(abs(gx), abs(gy), abs(gz))
        assert grad_max < Decimal("1e-90"), f"Stationary gradient too large: {grad_max}"

        lf = (
            dh * abs(x).ln()
            + dl * abs(Decimal(1) - x).ln()
            + dk * abs(y).ln()
            + ds * abs(Decimal(1) - y).ln()
            + dj * abs(z).ln()
            + dq * abs(Decimal(1) - z).ln()
            - dv13 * abs(D).ln()
        )
        if e == 0:
            c0 = -lf
            print(f"Inner (x0,y0,z0): grad_inf = {grad_max:.2e} < 1e-90 [PASS]")
        else:
            c1 = lf
            print(f"Outer (x1,y1,z1): grad_inf = {grad_max:.2e} < 1e-90 [PASS]")

    # 4. Build all 30 level-4 cosets and certify exact pair-sum splittings
    NUM, DEN = [], []
    for f in range(6):
        for u in range(FR[f]):
            NUM.append([theta(FN[f][i], u) for i in range(4)])
            DEN.append([theta(FD[f][i], u) for i in range(4)])
    assert len(NUM) == 30

    SP = ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))
    PR = []
    for c in range(30):
        found = None
        for n in SP:
            for d in SP:
                for sw in (0, 1):
                    d0, d1 = (2, 0) if sw else (0, 2)
                    n1 = V[NUM[c][n[0]]] + V[NUM[c][n[1]]]
                    n2 = V[NUM[c][n[2]]] + V[NUM[c][n[3]]]
                    e1 = V[DEN[c][d[d0]]] + V[DEN[c][d[d0 + 1]]]
                    e2 = V[DEN[c][d[d1]]] + V[DEN[c][d[d1 + 1]]]
                    if n1 == e1 and n2 == e2:
                        found = (
                            (NUM[c][n[0]], NUM[c][n[1]], DEN[c][d[d0]], DEN[c][d[d0 + 1]]),
                            (NUM[c][n[2]], NUM[c][n[3]], DEN[c][d[d1]], DEN[c][d[d1 + 1]]),
                        )
                        break
                if found:
                    break
            if found:
                break
        assert found is not None, f"No valid pair splitting for coset {c}"
        PR.append(found)

    aff = [set() for _ in range(16)]
    for c in range(30):
        for pair in PR[c]:
            for sy in pair:
                aff[sy].add(c)
    aff_list = [sorted(s) for s in aff]

    # 5. Exact integer breakpoint sweep over all distinct denominators dv
    dv_map = {}
    for i, val in enumerate(V):
        dv_map.setdefault(val, []).append(i)
    dv = sorted(dv_map.keys())
    who = [dv_map[val] for val in dv]
    nd = len(dv)
    a = [1] * nd
    F = [0] * 16
    cval = [0] * 30

    def recompute(c: int):
        p0, p1 = PR[c]
        v0 = F[p0[2]] + F[p0[3]] - F[p0[0]] - F[p0[1]]
        v1 = F[p1[2]] + F[p1[3]] - F[p1[0]] - F[p1[1]]
        assert -1 <= v0 <= 1 and -1 <= v1 <= 1, f"Lemma 4.1 violated at coset {c}: ({v0},{v1})"
        cval[c] = 2 if (v0 < 0 and v1 < 0) else (1 if (v0 + v1 < 0) else 0)

    for c in range(30):
        recompute(c)
    sv = max(cval)
    assert sv == 0, "Saving must start at 0 on [0, 1/M)"

    # Collect exact Fraction intervals [lo, hi) for Omega_E (sv >= 1) and Omega'_E (sv == 2)
    omega1_intervals = []
    omega2_intervals = []
    open1 = None
    open2 = None
    cut_points = 0
    changes = 0

    while True:
        best = -1
        ba, bd = 0, 1
        for z in range(nd):
            if a[z] >= dv[z]:
                continue
            if best < 0 or a[z] * bd < ba * dv[z]:
                best = z
                ba, bd = a[z], dv[z]
        if best < 0:
            break
        old = sv
        touched = set()
        for z in range(nd):
            if a[z] < dv[z] and a[z] * bd == ba * dv[z]:
                for sy in who[z]:
                    F[sy] += 1
                    touched.update(aff_list[sy])
                a[z] += 1
        for c in touched:
            recompute(c)
        sv = max(cval)
        cut_points += 1
        if sv != old:
            changes += 1
            om = Fraction(ba, bd)
            if old == 0 and sv >= 1:
                open1 = om
            elif old >= 1 and sv == 0:
                omega1_intervals.append((open1, om))
                open1 = None
            if old < 2 and sv == 2:
                open2 = om
            elif old == 2 and sv < 2:
                omega2_intervals.append((open2, om))
                open2 = None

    if open1 is not None:
        omega1_intervals.append((open1, Fraction(1, 1)))
    if open2 is not None:
        omega2_intervals.append((open2, Fraction(1, 1)))

    # 6. Certify Lemma 5.1 using EXACT rational arithmetic!
    first1 = omega1_intervals[0][0]
    first2 = omega2_intervals[0][0]
    assert first1 >= Fraction(1, M), f"Lemma 5.1 violated for Omega_E: {first1} < 1/{M}"
    assert first2 >= Fraction(1, N), f"Lemma 5.1 violated for Omega'_E: {first2} < 1/{N}"
    print(f"Lemma 5.1 EXACT: min(Omega_E) = {first1} >= 1/{M} [PASS]")
    print(f"Lemma 5.1 EXACT: min(Omega'_E) = {first2} >= 1/{N} [PASS]")
    print(f"Cut points = {cut_points}, changes = {changes}, components = {len(omega1_intervals)} / {len(omega2_intervals)}")

    # 7. Evaluate int_{Omega_E} dpsi and int_{Omega'_E} dpsi with rigorous error bounds
    # Cache psi values on exact Fractions so shared endpoints are computed once
    psi_cache = {}

    def get_psi(fr: Fraction):
        if fr not in psi_cache:
            x = Decimal(fr.numerator) / Decimal(fr.denominator)
            psi_cache[fr] = psi_dec(x)
        return psi_cache[fr]

    S1, err1 = Decimal(0), Decimal(0)
    for lo, hi in omega1_intervals:
        plo, elo = get_psi(lo)
        phi, ehi = get_psi(hi)
        S1 += phi - plo
        err1 += elo + ehi

    S2, err2 = Decimal(0), Decimal(0)
    for lo, hi in omega2_intervals:
        plo, elo = get_psi(lo)
        phi, ehi = get_psi(hi)
        S2 += phi - plo
        err2 += elo + ehi

    # Rigorous interval arithmetic for c2 and mu (include generous 1e-85 error on c0, c1)
    c_err = Decimal("1e-85")
    s_err = err1 + err2
    c2 = Decimal(M + N + Q) - (S1 + S2)
    assert c0 - c2 > Decimal(0), f"c0 <= c2 ({c0} <= {c2})"

    num_lo = (c0 + c1) - Decimal(2) * c_err
    num_hi = (c0 + c1) + Decimal(2) * c_err
    den_lo = (c0 - c2) - (c_err + s_err)
    den_hi = (c0 - c2) + (c_err + s_err)
    mu_lo = num_lo / den_hi
    mu_hi = num_hi / den_lo
    mu_mid = (c0 + c1) / (c0 - c2)

    dt = time.perf_counter() - t0
    print(f"c0                = {c0:.50f}")
    print(f"c1                = {c1:.50f}")
    print(f"int_{{Omega_E}}  dpsi = {S1:.50f}")
    print(f"int_{{Omega'_E}} dpsi = {S2:.50f}")
    print(f"c2                = {c2:.50f}")
    print(f"c0 - c2           = {c0 - c2:.50f}")
    print(f"Rigorous S error  < {s_err:.3e}")
    print(f"mu(zeta(3))       = {mu_mid:.60f}")
    print(f"Certified interval: [{mu_lo:.65f},")
    print(f"                     {mu_hi:.65f}]")
    print(f"Interval width    = {mu_hi - mu_lo:.3e}  (elapsed {dt:.2f}s)")
    return mu_mid, mu_lo, mu_hi


if __name__ == "__main__":
    if len(sys.argv) == 9:
        pts = tuple(int(x) for x in sys.argv[1:9])
        certify(pts, "CLI CANDIDATE")
    else:
        certify((16, 17, 19, 15, 12, 11, 9, 13), "SELF-TEST: Rhin-Viola (2001)")
        certify((1167, 1240, 1386, 1094, 875, 802, 656, 948), "NEW CANDIDATE: 73*RV - 1")
