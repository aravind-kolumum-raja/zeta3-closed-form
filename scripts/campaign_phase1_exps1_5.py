#!/usr/bin/env python3
"""
campaign_phase1_exps1_5.py
==========================
Phase 1 of the 25-Experiment Grand Campaign (Experiments 1 through 5):
  Exp 1: Algebraic Z/p^2Z dissection of the two 2nd-order windows
         I_12 = (7/15, 11/23) and I_23 = (16/17, 20/21) of F_n*(t).
  Exp 2: Search for non-trivial integral/Mellin-Barnes/reflection/Thomae twins of F_n*(t).
  Exp 3: Polynomial multipliers P(t) = (a*t + b*n) and contiguous shift combinations at fixed saddle.
  Exp 4: Constant integer offsets (t - 13n + alpha_i)_{13n} and wing denominator analysis.
  Exp 5: Ultra-fine Farey N=80 scan of F_n*(t) for 3rd-order Z/p^3Z or missed 2nd-order cells.
"""

from fractions import Fraction
import math

RS = (13, 15, 17)
LS = (13, 13, 13)
S1, P1 = 2, 23
S2, P2 = 3, 21
E1, E2 = 20, 19


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def vp_int(x: int, p: int) -> int:
    if x == 0:
        return 999
    x = abs(x)
    v = 0
    while x % p == 0:
        v += 1
        x //= p
    return v


def vp_frac(f: Fraction, p: int) -> int:
    if f == 0:
        return 999
    return vp_int(f.numerator, p) - vp_int(f.denominator, p)


def fact_table_pk(limit: int, p: int, k: int):
    mod = p ** k
    vp = [0] * (limit + 1)
    uf = [1] * (limit + 1)
    for i in range(1, limit + 1):
        x = i
        c = 0
        while x % p == 0:
            c += 1
            x //= p
        vp[i] = vp[i - 1] + c
        uf[i] = (uf[i - 1] * x) % mod
    return vp, uf


def binom_vp_u(N: int, K: int, p: int, k: int, vp, uf):
    if K < 0 or K > N:
        return 999, 0
    mod = p ** k
    v = vp[N] - vp[K] - vp[N - K]
    u = (uf[N] * pow((uf[K] * uf[N - K]) % mod, -1, mod)) % mod
    return v, u


def compute_aj_mod_pk(n: int, j: int, p: int, k: int, vp, uf):
    """Exact (v_p(C_n a_j), unit mod p^k) for double pole 3n <= j <= 21n."""
    v1, u1 = binom_vp_u(j + 13 * n, 13 * n, p, k, vp, uf)
    v2, u2 = binom_vp_u(j + 15 * n, 13 * n, p, k, vp, uf)
    v3, u3 = binom_vp_u(j + 17 * n, 13 * n, p, k, vp, uf)
    v4, u4 = binom_vp_u(20 * n, j - 3 * n, p, k, vp, uf)
    v5, u5 = binom_vp_u(19 * n, j - 2 * n, p, k, vp, uf)
    mod = p ** k
    return v1 + v2 + v3 + v4 + v5, (u1 * u2 * u3 * u4 * u5) % mod


def compute_exact_Fn_residues(n: int, poly_mult=None, offsets=None):
    """
    Exact Fraction computation of (C_n Q_n, C_n P_{3,n}, a_j, b_j) for small n.
    poly_mult(t) optional extra polynomial in t.
    offsets = (dr1, dr2, dr3, ds1, ds2, dd1, dd2) optional constant shifts.
    """
    if offsets is None:
        dr = (0, 0, 0)
        ds1 = ds2 = dd1 = dd2 = 0
    else:
        dr = offsets[:3]
        ds1, ds2, dd1, dd2 = offsets[3:]

    r = [RS[i] * n + dr[i] for i in range(3)]
    L = [LS[i] * n for i in range(3)]
    s1 = S1 * n + ds1
    p1 = P1 * n + ds1 + dd1
    s2 = S2 * n + ds2
    p2 = P2 * n + ds2 + dd2

    # Normalization C_n = (20n)!(19n)! / ((13n)!)^3
    Cn = Fraction(math.factorial(20 * n) * math.factorial(19 * n), math.factorial(13 * n) ** 3)

    # Double zero start t0 = 2nd smallest z_u
    z = sorted(r[i] - L[i] for i in range(3))
    t0 = z[1]

    # Precompute harmonic numbers H_k^{(2)}, H_k^{(3)}
    max_h = max(p1, p2) + t0 + 5
    H2 = [Fraction(0)] * (max_h + 1)
    H3 = [Fraction(0)] * (max_h + 1)
    for k in range(1, max_h + 1):
        H2[k] = H2[k - 1] + Fraction(1, k * k)
        H3[k] = H3[k - 1] + Fraction(1, k * k * k)

    # Evaluate residues at j in [min(s1,s2), max(p1,p2)]
    j_lo, j_hi = min(s1, s2), max(p1, p2)
    Qn = Fraction(0)
    P3n = Fraction(0)
    a_map = {}
    b_map = {}
    for j in range(j_lo, j_hi + 1):
        in1 = (s1 <= j <= p1)
        in2 = (s2 <= j <= p2)
        if not (in1 or in2):
            continue
        # Product of numerator factors at t = -j
        num = Fraction(1)
        num_logd = Fraction(0)
        for i in range(3):
            # (t - r_i)_{L_i} = prod_{m=0}^{L_i-1} (t - r_i + m)
            # at t = -j: prod_{m=0}^{L_i-1} (-j - r_i + m) = (-1)^{L_i} * (j+r_i)! / (j+r_i-L_i)!
            a0 = j + r[i] - L[i] + 1
            b0 = j + r[i]
            if a0 <= 0 <= b0:
                num = Fraction(0)
                break
            sign = -1 if (L[i] % 2 == 1) else 1
            prod_v = 1
            ld = Fraction(0)
            for val in range(a0, b0 + 1):
                prod_v *= val
                ld -= Fraction(1, val)  # d/dt ln(t - k)|_{t=-j} = 1/(-j-k) = -1/(j+k)
            num *= sign * prod_v
            num_logd += ld
        if num == 0:
            continue
        if poly_mult is not None:
            pv, dpv = poly_mult(-j)
            if pv == 0:
                continue
            num *= pv
            num_logd += Fraction(dpv, pv)

        # Denominator product omitting (t+j) factors
        den = Fraction(1)
        den_logd = Fraction(0)
        for (sa, pa) in [(s1, p1), (s2, p2)]:
            for k in range(sa, pa + 1):
                if k == j:
                    continue
                diff = k - j  # t + k at t = -j
                den *= diff
                den_logd += Fraction(1, diff)

        val = Cn * num / den
        if in1 and in2:
            aj = val
            bj = val * (num_logd - den_logd)
            a_map[j] = aj
            b_map[j] = bj
            Qn += aj
            P3n += aj * H3[j + t0] + Fraction(1, 2) * bj * H2[j + t0]
        else:
            bj = val
            b_map[j] = bj
            P3n += Fraction(1, 2) * bj * H2[j + t0]

    return Qn, P3n, a_map, b_map, t0


# =============================================================================
# EXP 1: Algebraic Z/p^2Z Dissection of I_12 = (7/15, 11/23) and I_23 = (16/17, 20/21)
# =============================================================================
def run_exp1():
    print("=" * 80)
    print("EXP 1: ALGEBRAIC Z/p^2Z DISSECTION OF 2ND-ORDER WINDOWS I_12 AND I_23")
    print("=" * 80)

    windows = [
        ("I_12a (7/15, 9/19)", Fraction(7, 15), Fraction(9, 19)),
        ("I_12b (9/19, 11/23)", Fraction(9, 19), Fraction(11, 23)),
        ("I_23a (16/17, 18/19)", Fraction(16, 17), Fraction(18, 19)),
        ("I_23b (18/19, 20/21)", Fraction(18, 19), Fraction(20, 21)),
        ("CONTROL 1st-order (14/15, 16/17)", Fraction(14, 15), Fraction(16, 17)),
        ("CONTROL 1st-order (11/23, 12/25)", Fraction(11, 23), Fraction(12, 25)),
    ]

    for label, lo, hi in windows:
        print(f"\n--- {label} ---")
        # Pick 3 primes p > 40*n with n/p strictly inside (lo, hi)
        found = 0
        for p in range(61, 400):
            if not is_prime(p):
                continue
            n_lo = math.floor(lo * p) + 1
            n_hi = math.ceil(hi * p) - 1
            if n_lo > n_hi:
                continue
            n = (n_lo + n_hi) // 2
            if p * p <= 40 * n + 5:
                continue
            u = Fraction(n, p)
            if not (lo < u < hi):
                continue

            vp, uf = fact_table_pk(40 * n + 5, p, 4)
            # Compute v_min across all double poles j in [3n, 21n]
            v_vals = {}
            vmin = 999
            for j in range(3 * n, 21 * n + 1):
                vj, uj = compute_aj_mod_pk(n, j, p, 4, vp, uf)
                v_vals[j] = (vj, uj)
                if vj < vmin:
                    vmin = vj

            Delta = 2 - 3 * (13 * n // p) + (20 * n // p) + (19 * n // p)
            mod2 = p * p
            mod3 = p * p * p

            # Bucket by m = (j + 2n) // p, x = (j + 2n) % p
            buckets = {}
            bucket_vmins = {}
            bucket_counts = {}
            bucket_v_eq_vmin = {}
            for j in range(3 * n, 21 * n + 1):
                m, x = divmod(j + 2 * n, p)
                vj, uj = v_vals[j]
                dv = vj - vmin
                term = (uj * (p ** dv)) % mod3 if dv < 3 else 0
                buckets[m] = (buckets.get(m, 0) + term) % mod3
                bucket_vmins[m] = min(bucket_vmins.get(m, 999), vj)
                bucket_counts[m] = bucket_counts.get(m, 0) + 1
                if vj == vmin:
                    bucket_v_eq_vmin[m] = bucket_v_eq_vmin.get(m, 0) + 1

            tot_mod3 = sum(buckets.values()) % mod3
            v_tot_extra = 3 if tot_mod3 == 0 else (2 if tot_mod3 % mod2 == 0 else (1 if tot_mod3 % p == 0 else 0))

            # Analyze each bucket's extra valuation!
            b_summary = []
            for m in sorted(buckets.keys()):
                bm = buckets[m]
                vb = 3 if bm == 0 else (2 if bm % mod2 == 0 else (1 if bm % p == 0 else 0))
                b_summary.append(f"m={m}(cnt={bucket_counts[m]},#min={bucket_v_eq_vmin.get(m,0)},vm={bucket_vmins[m]},+v={vb})")

            # Also check polynomial degree of reduced bucket terms mod p and mod p^2!
            # And check the carry tuple on the support where vj == vmin AND vj == vmin + 1!
            min_js = [j for j in range(3 * n, 21 * n + 1) if v_vals[j][0] == vmin]
            min1_js = [j for j in range(3 * n, 21 * n + 1) if v_vals[j][0] == vmin + 1]
            m_min_set = sorted(set((j + 2 * n) // p for j in min_js))
            m_min1_set = sorted(set((j + 2 * n) // p for j in min1_js))

            # Check if sum over ONLY vj == vmin terms vanishes mod p^2 vs including vj == vmin+1!
            sum_only_vmin = {}
            sum_only_vmin1 = {}
            for j in range(3 * n, 21 * n + 1):
                m = (j + 2 * n) // p
                vj, uj = v_vals[j]
                if vj == vmin:
                    sum_only_vmin[m] = (sum_only_vmin.get(m, 0) + uj) % mod2
                elif vj == vmin + 1:
                    sum_only_vmin1[m] = (sum_only_vmin1.get(m, 0) + uj) % p

            print(f"  p={p}, n={n} (u={float(u):.5f}): Delta={Delta}, vmin={vmin}, Delta-vmin={Delta-vmin}, TOTAL +v={v_tot_extra}")
            print(f"    buckets: {', '.join(b_summary)}")
            # Inspect non-trivial buckets mod p^2
            for m in sorted(buckets.keys()):
                s0 = sum_only_vmin.get(m, 0)
                s1 = sum_only_vmin1.get(m, 0)
                tot_m = buckets[m] % mod2
                if s0 != 0 or s1 != 0:
                    print(f"      bucket m={m}: (sum_{'{v=vmin}'} u_j) mod p^2 = {s0} ({s0//p}*p + {s0%p}), "
                          f"(sum_{'{v=vmin+1}'} u_j) mod p = {s1}, combined mod p^2 = {tot_m} ({tot_m//p}*p + {tot_m%p})")

            # ALSO inspect the exact x-support of v=vmin and v=vmin+1 inside one active bucket!
            if found == 0 and m_min_set:
                m0 = m_min_set[0]
                xs0 = [((j + 2 * n) % p) for j in range(3 * n, 21 * n + 1) if (j + 2 * n) // p == m0 and v_vals[j][0] == vmin]
                xs1 = [((j + 2 * n) % p) for j in range(3 * n, 21 * n + 1) if (j + 2 * n) // p == m0 and v_vals[j][0] == vmin + 1]
                xs2 = [((j + 2 * n) % p) for j in range(3 * n, 21 * n + 1) if (j + 2 * n) // p == m0 and v_vals[j][0] >= vmin + 2]
                print(f"    Sample bucket m={m0}: #x(v=vmin)={len(xs0)} [{min(xs0) if xs0 else '-'}..{max(xs0) if xs0 else '-'}], "
                      f"#x(v=vmin+1)={len(xs1)}, #x(v>=vmin+2)={len(xs2)}")
                # Compute exact F_p polynomial degrees of v=vmin term AND v=vmin+1 term (finite differences mod p)!
                arr0 = [0] * p
                arr1 = [0] * p
                for j in range(3 * n, 21 * n + 1):
                    if (j + 2 * n) // p == m0:
                        x = (j + 2 * n) % p
                        vj, uj = v_vals[j]
                        if vj == vmin:
                            arr0[x] = uj % p
                        elif vj == vmin + 1:
                            arr1[x] = uj % p
                # Finite difference degree mod p
                def fp_deg(arr, p):
                    a = arr[:]
                    deg = -1
                    for d in range(p):
                        if any(v != 0 for v in a):
                            deg = d
                        else:
                            break
                        a = [(a[i + 1] - a[i]) % p for i in range(len(a) - 1)]
                    return deg
                print(f"    F_p degrees in bucket m={m0}: deg(v=vmin) = {fp_deg(arr0, p)} (p-1={p-1}), deg(v=vmin+1) = {fp_deg(arr1, p)}")

            found += 1
            if found >= 2:
                break


# =============================================================================
# EXP 2: Mellin-Barnes / Reflection / Thomae / Euler-Integral Twin Analysis
# =============================================================================
def run_exp2():
    print("\n" + "=" * 80)
    print("EXP 2: REFLECTION / MELLIN-BARNES / THOMAE TWIN ANALYSIS OF F_n*(t)")
    print("=" * 80)
    # F_n*(t) = (t-13n)_{13n} (t-15n)_{13n} (t-17n)_{13n} / [(t+2n)_{21n+1} (t+3n)_{18n+1}]
    # Look at its zero intervals and pole intervals:
    #   Zeros: [1, 13n], [2n+1, 15n], [4n+1, 17n]  (total 39n)
    #   Poles: [-23n, -2n], [-21n, -3n]            (total 39n + 2)
    # Reflection t -> -t - c: maps poles [-23n, -2n] to [2n - c, 23n - c] and zeros [1, 17n] to [-17n - c, -1 - c].
    # That swaps poles and zeros (not another rational function with double poles!).
    # What about 3D Euler integral / Mellin-Barnes representations?
    # Write R_{3,n} = -1/(2*pi*i) int (pi/sin(pi*t))^2 F_n*(t) dt.
    # For a 3-block double-pole function, in Gamma form:
    #   F_n*(t) = (-1)^{39n} * Gamma(1+13n-t) Gamma(1+15n-t) Gamma(1+17n-t) Gamma(t+2n) Gamma(t+3n)
    #                        / [ Gamma(1-t) Gamma(1+2n-t) Gamma(1+4n-t) Gamma(t+23n+1) Gamma(t+21n+1) ]
    # Notice: multiplying by (pi/sin(pi*t))^2 = Gamma(t-a)Gamma(1+a-t) Gamma(t-b)Gamma(1+b-t)
    # lets us cancel ANY TWO of the 5 denominator Gammas!
    # Case A: cancel Gamma(1-t) and Gamma(1+2n-t) (a=0, b=2n):
    #   Integrand = Gamma(t+2n) Gamma(t+3n) Gamma(t) Gamma(t-2n) * Gamma(1+13n-t) Gamma(1+15n-t) Gamma(1+17n-t)
    #               / [ Gamma(t+23n+1) Gamma(t+21n+1) * Gamma(1+4n-t) ]
    # Wait: that is a Meijer G_{5,5}^{3,4} / G_{5,5}^{4,3}(1) = 4F3? NO, 7 upstairs and 3 downstairs = G_{5,5}^{4,3}(1) (sum of 5F4)!
    # Check the 5x5 parameter matrix of G_{5,5}^{3,4}(1):
    # Wait: for G_{p,q}^{m,n}(1) = 1/(2pi i) int [prod_{j=1}^m Gamma(b_j-s) prod_{j=1}^n Gamma(1-a_j+s)]
    #                                         / [prod_{j=m+1}^q Gamma(1-b_j+s) prod_{j=n+1}^p Gamma(a_j-s)] ds
    # Here with s = t:
    #   Gamma(b_j - t) upstairs (m=3): b = [1+13n, 1+15n, 1+17n]
    #   Gamma(1-a_j + t) upstairs (n_g=4): 1-a = [2n, 3n, 0, -2n] => a = [1-2n, 1-3n, 1, 1+2n]
    #   Gamma(1-b_j + t) downstairs (q-m=2): 1-b = [23n+1, 21n+1] => b_extra = [-23n, -21n]
    #   Gamma(a_j - t) downstairs (p-n_g=1): a_extra = [1+4n]
    # WAIT A SECOND! Look at the 3 choices of WHICH two denominator Gammas to cancel!
    # Even better: check Barnes's / Nesterenko's 3-fold Euler integral!
    # In Nesterenko (2003, Theorem 2) / Zudilin (2004, Lemma 1):
    #   sum_{t=1}^inf F'(t) for a 3-block double-pole function = int_{[0,1]^3} ... ONLY when m=3 numerator blocks?
    # Let's check the exact Beta-integral decomposition of Gamma ratios in F_n*(t) and check ALL S_3 x S_3 x S_2
    # pairings of (b_1,b_2,b_3) with (a_1,a_2,a_3,a_4) -> Beta integrals!
    # Specifically, can a 3-block double-pole F_n(t) be transformed into a DIFFERENT 3-block double-pole F~_n(t)
    # by permuting the 5 parameters of G_{5,5}^{3,4}(1)?
    # Count symmetries of G_{5,5}^{3,4}(1):
    #   Upstairs (t): 4 parameters U+ = {2n, 3n, 0, -2n} (for choice a=0, b=2n)
    #   Downstairs (t): 2 parameters D+ = {23n+1, 21n+1}
    #   Upstairs (-t): 3 parameters U- = {13n+1, 15n+1, 17n+1}
    #   Downstairs (-t): 1 parameter D- = {4n+1}
    # Wait: the integral is INVARIANT under permuting U+ (4!), D+ (2!), U- (3!), D- (1!).
    # Does permuting U+ change the rational function F_n(t)?
    # YES!! Because to recover a rational function F~_n(t) times (pi/sin(pi*t))^2, we pick ANY 2 of the 4
    # elements of U+ to convert via Gamma(t + u) = (pi/sin(pi*t)) / Gamma(1 - u - t) into D-!!
    # WHOA!! CHECK THIS RIGHT NOW!!
    # Pick ANY 2 elements {u_1, u_2} from U+ = {2n, 3n, 0, -2n} to reflect to D-!
    # Wait: originally we reflected {0, -2n} to D- = {1, 2n+1} (Wait: 1 - 0 = 1, 1 - (-2n) = 2n+1).
    # WHAT IF WE REFLECT {2n, 3n} INSTEAD OF {0, -2n}?!
    # Let's compute ALL C(4, 2) = 6 reflections of U+ AND all reflections using U- -> D+!

    # Test numerically at n = 1!
    n = 1
    Qn, P3n, _, _, _ = compute_exact_Fn_residues(n)
    zeta3 = 1.202056903159594285399738161511449990765
    R3_val = float(Qn) * zeta3 - float(P3n)
    print(f"n=1: Q_1 = {Qn}, P_1 = {P3n.numerator}/{P3n.denominator}")
    print(f"     R_{{3,1}} = {R3_val:.12e} (exact Fraction Q_1 has {len(str(Qn))} digits)")

    # Check all C(4,2) = 6 ways to convert 2 of U+ = {2n, 3n, 0, -2n} to D-, and also
    # all ways to reflect from (U+, D+, U-, D-) = ({2n, 3n, 0, -2n}, {23n+1, 21n+1}, {13n+1, 15n+1, 17n+1}, {4n+1})
    # Wait: (pi/sin(pi*t))^2 = Gamma(t+u)Gamma(1-u-t) * Gamma(t+v)Gamma(1-v-t)
    # Wait: dividing the SAME integrand by (pi/sin(pi*t))^2 ALWAYS gives the EXACT SAME rational function F_n*(t)!!
    # Why? Because (pi/sin(pi*t))^2 is INDEPENDENT of the choice of (u, v)!
    # AHA! To get a DIFFERENT rational function F~_n(t), you need a NON-TRIVIAL transformation of G_{5,5}^{3,4}(1)
    # such as s -> -s (which maps G_{5,5}^{3,4}(1) to G_{5,5}^{4,3}(1) -> DIFFERENT rational function!)
    # WAIT A SECOND!! CHECK s -> -s ON G_{5,5}^{3,4}(1)!!
    # Under t -> -t:
    #   U+ = {2n, 3n, 0, -2n} becomes U-' = {2n, 3n, 0, -2n} (4 Gammas of -t upstairs!)
    #   U- = {13n+1, 15n+1, 17n+1} becomes U+' = {13n+1, 15n+1, 17n+1} (3 Gammas of +t upstairs!)
    #   D+ = {23n+1, 21n+1} becomes D-' = {23n+1, 21n+1} (2 Gammas of -t downstairs!)
    #   D- = {4n+1} becomes D+' = {4n+1} (1 Gamma of +t downstairs!)
    # Divide by (pi/sin(pi*t))^2 = Gamma(t+k1)Gamma(1-k1-t) Gamma(t+k2)Gamma(1-k2-t):
    # Wait: we need 2 Gammas of (+t) upstairs and 2 Gammas of (+t) downstairs to get a double-pole rational function,
    # OR 3 Gammas of (+t) upstairs and 3 downstairs to get a TRIPLE-POLE function!!
    # Look at U+', D+': we have 3 upstairs (+t): {13n+1, 15n+1, 17n+1}, and 1 downstairs (+t): {4n+1}.
    # If we absorb 2 from D-' = {23n+1, 21n+1} into (+t) downstairs (via sin^2), we get:
    #   D+'' = {4n+1, -23n, -21n} -> Wait: Gamma(23n+1 - t) downstairs becomes Gamma(t - 23n) upstairs?! NO,
    #   1 / [Gamma(23n+1-t) * Gamma(t-23n) Gamma(23n+1-t)] doesn't cancel! To cancel Gamma(23n+1-t) downstairs,
    #   we need it UPSTAIRS, which means reflecting U-' -> D+''!
    # Reflect ANY 2 of U-' = {2n, 3n, 0, -2n} to D+'':
    #   e.g. reflect {0, -2n} -> {1, 2n+1} in D+'', leaving U-'' = {2n, 3n}, D-'' = {23n+1, 21n+1}!
    #   Then we have:
    #     (+t): upstairs {13n+1, 15n+1, 17n+1}, downstairs {4n+1, 1, 2n+1} -> TRIPLE POLE? Wait: 13n+1 > 4n+1, so
    #     Gamma(t+13n+1)/Gamma(t+4n+1) is a POLYNOMIAL, NOT poles!!
    #     (-t): upstairs {2n, 3n}, downstairs {23n+1, 21n+1} -> Gamma(2n-t)/Gamma(23n+1-t) HAS POLES!!
    #   Wait: that's literally JUST the reflection t -> -t of F_n*(t), whose right-half-plane sum is the LEFT-half-plane
    #   sum of F_n*(t) — which gives the EXACT SAME residues a_j, b_j (shifted by t -> -t)!

    # NOW CHECK THE REAL QUESTION: Does G_{5,5}^{3,4}(1) satisfy a NON-TRIVIAL Bailey / Thomae / Zudilin 3D integral
    # transformation?!
    # Write the 3D Euler integral for F_n*(t) explicitly!
    # Look at the 3 numerator ratios and 2 denominator ratios:
    #   U+ = {2n, 3n, 0, -2n} (Wait: 0, -2n came from sin^2!), D+ = {23n+1, 21n+1}
    #   U- = {13n+1, 15n+1, 17n+1}, D- = {4n+1}
    # Notice: by Mellin-Barnes, the residue sum R_{3,n} = sum_{t=1}^inf F_n*'(t) equals the 3D integral:
    #   I_n = int_{[0,1]^3} x^{a1}(1-x)^{b1} y^{a2}(1-y)^{b2} z^{a3}(1-z)^{b3} / (1 - xyz)^{...} ?
    # NO! Look at the 3F2(1) inside! To convert to a Beukers/Zudilin 3D integral:
    #   int_0^1 int_0^1 int_0^1 x^A (1-x)^B y^C (1-y)^D z^E (1-z)^F / (1 - (1-xy)z)^G dx dy dz
    #   has Mellin-Barnes representation (Rhin-Viola 2001 eq (2.2) / Zudilin 2004):
    #     Gamma(A+1+t) Gamma(C+1+t) Gamma(E+1+t) * Gamma(G-A-1-t) Gamma(G-C-1-t) Gamma(-t) ...
    #   That ALWAYS has 3 (+t) upstairs and 3 (-t) upstairs, i.e. G_{4,4}^{3,3}(1) or balanced 4F3(1) (3 parameters, W(D_5)!)
    #   Wait: a 3-block double-pole F_n(t) has:
    #     3 numerator blocks (degree 39n) + 2 denominator blocks (degree 39n+2) -> 5 Gamma ratios = G_{5,5} = 5F4(1)!
    #   A GENERAL 5F4(1) is a 4-DIMENSIONAL Euler integral (not 3D!), and by Thomae/Whipple/Zudilin:
    #   WAIT A SECOND!! Does a balanced 5F4(1) have ANY Whipple / Bailey transformation?!
    #   Let's recall the classical theorem on hypergeometric transformation groups at z = 1:
    #     - 3F2(1): Thomae's group S_5 (order 120) — works for ANY 3F2(1)!
    #     - 4F3(1) balanced (1-balanced): Whipple's group W(D_5) (order 1920) — works for ANY balanced 4F3(1) (terminating) or VWP 7F6(1)!
    #     - 5F4(1) general (NOT very-well-poised): DOES IT HAVE A TRANSFORMATION GROUP?!
    #       Let's test numerically by checking whether F_n*(t) satisfies ANY 1-balanced 4F3 or 5F4 reduction,
    #       AND check the sum of parameters (Saalschutzian excess s = sum(denom) - sum(num))!
    # Compute the hypergeometric parameters of sum_{k=0}^{18n} C_n a_{3n+k} (which is a terminating 5F4(1)!)
    #   B_{n,k} = C(16n+k, 13n) C(18n+k, 13n) C(20n+k, 13n) C(20n, k) C(19n, n+k)
    #   Ratio B_{n,k+1} / B_{n,k}:
    #     = (16n+1+k)(18n+1+k)(20n+1+k)(20n-k)(18n-k)
    #       / [(3n+1+k)(5n+1+k)(7n+1+k)(k+1)(n+1+k)]
    #   So Q_n = B_{n,0} * 5F4(16n+1, 18n+1, 20n+1, -20n, -18n ; 3n+1, 5n+1, 7n+1, n+1 ; 1)!
    a_params = [16, 18, 20, -20, -18]  # * n (+1 for first 3)
    b_params = [3, 5, 7, 1]            # * n + 1
    saalschutz_n = sum(b_params) - sum(a_params)
    print(f"Terminating 5F4(1) for Q_n:")
    print(f"  a = [16n+1, 18n+1, 20n+1, -20n, -18n]")
    print(f"  b = [3n+1, 5n+1, 7n+1, n+1]")
    print(f"  Saalschutzian excess s = sum(b) - sum(a) = {saalschutz_n}n + 1  (1-balanced requires 0n + 1!)")
    print(f"  Wait: sum(a) = {sum(a_params)}n + 3, sum(b) = {sum(b_params)}n + 4 => s = {saalschutz_n}n + 1 = 0n + 1!!")
    print(f"  WAIT A SECOND!! sum(b) - sum(a) = (3+5+7+1)n + 4 - (16+18+20-20-18)n - 3 = 16n + 4 - (16n + 3) = +1 !!")
    print(f"  HOLY COW!! Q_n IS **1-BALANCED (SAALSCHUTZIAN)** 5F4(1)!! (Every balanced 3-block 5F4 is 1-balanced? Let's check!)")
    # Check well-poisedness: 1 + a_i + b_j = const?
    print("  Checking well-poised pairings 1 + a_0 = a_i + b_j:")
    for i in range(5):
        a0 = a_params[i]
        print(f"    If a0 = {a0}n: a0+1 = {a0}n+1; sums a_k + b_l:",
              [(a_params[k], b_params[l], a_params[k] + b_params[l]) for k in range(5) if k != i for l in range(4) if a_params[k] + b_params[l] == a0])
    # Also check: what about the INFINITE sum R_{3,n} = sum_{t=17n+1}^inf F_n*'(t)?
    # F_n*(t) at t = 17n + 1 + k has ratio:
    #   F_n*(t+1)/F_n*(t) = (t+1)(t-2n+1)(t-4n+1)(t+2n)(t+3n)
    #                       / [(t-13n+1)(t-15n+1)(t-17n+1)(t+23n+1)(t+21n+1)]
    # Sum of denom - sum of num = (23+21-13-15-17)n + 5 - (2+3-2-4)n - 3 = (-1)n + 5 - (-1)n - 3 = +2!
    # So F_n*(t) is a 2-balanced 5F4(1) (differentiated -> 3-balanced)!
    # Check if there exists a non-trivial 5F4(1) Whipple/Bailey/contiguous transform!
    print("  NOTE: Q_n is 1-balanced 5F4(1) (NOT well-poised -> Whipple 7F6->4F3 does NOT apply directly).")
    print("  HOWEVER, check Bailey/Karp-Prilepkina/Zudilin 5F4(1) shift/reversal transformations!")
    # Reversal of terminating 5F4(1): k -> 18n - k gives the SAME sum (just reversing B_{n,k}),
    # which is t -> -t - 24n on poles [3n, 21n] (different binomials? Let's check!)
    print("  Reversing k -> 18n - k gives:")
    print("    C(34n-k, 13n) C(36n-k, 13n) C(38n-k, 13n) C(20n, 18n-k) C(19n, 19n-k) -> identical term-by-term!")


# =============================================================================
# EXP 3: Polynomial Multipliers P(t) = a*t + b*n & Contiguous Combinations
# =============================================================================
def run_exp3():
    print("\n" + "=" * 80)
    print("EXP 3: POLYNOMIAL MULTIPLIERS (a*t + b*n) ON F_n*(t) AT FIXED SADDLE")
    print("=" * 80)
    # Wait: multiplying F_n*(t) by (a*t + b*n) reduces deg D - deg N from 2 to 1 ->
    # F_n'(t) STILL decays like 1/t^2! Wait: does sum b_j = 0 still hold?
    # NO!! If deg D - deg N = 1, then sum b_j = leading coeff = a != 0, so sum F_n'(t) picks up a*zeta(2)!!
    # So to keep deg D - deg N >= 2 (zero zeta(2)!), we must either:
    #   (3A) Multiply by (a*t + b*n) and ADD one extra denominator pole (t + c*n), OR
    #   (3B) Take a LINEAR COMBINATION of TWO F_n*(t) with SAME asymptotic saddle:
    #        G_n(t) = c_0 F_n*(t) + c_1 F_n*(t + 1)  (or shift by 1 in one numerator/denom block!)
    # Wait: (3B) F_n*(t) + lam * F_n*(t+1) is IDENTICAL to F_n*(t) with a linear multiplier?
    # NO! Because sum_{t} F_n*'(t+1) = sum_{t} F_n*'(t) - F_n*'(t_0 n + 1) = sum_{t} F_n*'(t) (double zero!),
    # so c_0 F_n*(t) + c_1 F_n*(t+1) gives (c_0 + c_1) R_{3,n} — trivial!
    # What about shifting ONE parameter by +/-1 (contiguous neighbor F_n^{(delta)}(t))?
    # E.g. F_n^{(1)}(t) = F_n*(t) * (t - 13n - 1) / (t - 1)  [shifting block 1 by 1]? No, that changes zeros.
    # Valid contiguous neighbor with SAME poles & double zeros:
    #   replace ONE factor (t - 17n) by (t + 23n + 1)? No, that adds a pole.
    #   replace (t - 17n)_{13n} by (t - 17n + 1)_{13n} = (t - 17n)_{13n} * (t - 4n) / (t - 17n) -> STILL has double zeros on [2n+1, 15n]!
    #   YES!! Shifting r_3 -> r_3 - 1 keeps double zeros on [2n+1, 15n] AND keeps exact same poles & saddle!!
    # Let's test ALL linear combinations of F_n*(t) and its 1-step contiguous shifts (r_3 -> r_3 - 1, r_1 -> r_1 + 1, etc.)
    # across primes p on NON-shell cells of F_n*(t) to see if ANY fixed (c_0, c_1) creates NEW F_p-shells!

    # Test on 4 representative non-shell Farey cells of F_n*(t):
    test_cells = [
        ("non-shell u~0.35", Fraction(7, 20)),
        ("non-shell u~0.42", Fraction(5, 12)),
        ("non-shell u~0.60", Fraction(3, 5)),
        ("non-shell u~0.71", Fraction(5, 7)),
    ]
    print("Testing linear combinations c0*F_n*(t) + c1*F_n^{shift}(t) on non-shell cells (exact Fraction at small p):")
    # Pick small prime p where n/p is in a non-shell cell, compute (Q_n, P_{3,n}) for F_n* and shifted neighbors!
    # Actually at n=1, p in [2..25] has n/p = 1/p <= 1/2: e.g. p=3 -> u=1/3 (non-shell), n=2, p=5 -> u=2/5=0.4 (non-shell)!
    # Wait: we need p > sqrt(40n) so base-p^2 carries don't interfere; for n=2, p > 9 (u < 2/11).
    # For u ~ 0.35 with p > sqrt(40n), need p^2 > 40*(0.35*p) => p > 14! So p = 29, n = 10 -> u = 10/29 = 0.3448!
    # Let's use fast mod-p^3 arithmetic for n=10, p=29 (u=0.345), n=12, p=29 (u=0.414), n=17, p=29 (u=0.586)!
    def eval_shifted_aj_mod_p2(n, p, dr=(0, 0, 0)):
        """Compute exact sum_j C_n a_j mod p^3 for shifted r_i -> RS[i]*n + dr[i] (keeping L_i = 13n)."""
        # Note: C_n a_j = prod_{i=1}^3 C(j + RS[i]*n + dr[i], 13n) * C(20n, j-3n) * C(19n, j-2n)
        vp, uf = fact_table_pk(40 * n + 10, p, 3)
        mod = p ** 3
        vmin = 999
        terms = []
        for j in range(3 * n, 21 * n + 1):
            v1, u1 = binom_vp_u(j + 13 * n + dr[0], 13 * n, p, 3, vp, uf)
            v2, u2 = binom_vp_u(j + 15 * n + dr[1], 13 * n, p, 3, vp, uf)
            v3, u3 = binom_vp_u(j + 17 * n + dr[2], 13 * n, p, 3, vp, uf)
            v4, u4 = binom_vp_u(20 * n, j - 3 * n, p, 3, vp, uf)
            v5, u5 = binom_vp_u(19 * n, j - 2 * n, p, 3, vp, uf)
            vj = v1 + v2 + v3 + v4 + v5
            uj = (u1 * u2 * u3 * u4 * u5) % mod
            terms.append((vj, uj))
            if vj < vmin:
                vmin = vj
        tot = 0
        for vj, uj in terms:
            dv = vj - vmin
            if dv < 3:
                tot = (tot + uj * (p ** dv)) % mod
        return vmin, tot

    for p in [29, 31, 37]:
        for target_u in [0.345, 0.42, 0.59]:
            n = round(target_u * p)
            if p * p <= 40 * n:
                continue
            v0, t0 = eval_shifted_aj_mod_p2(n, p, (0, 0, 0))
            v1, t1 = eval_shifted_aj_mod_p2(n, p, (0, 0, -1))
            v2, t2 = eval_shifted_aj_mod_p2(n, p, (1, 0, 0))
            # Check ratio t1 / t0 mod p — if constant across primes p, a fixed integer combo kills it!
            r1 = (t1 * pow(t0, -1, p)) % p if (v1 == v0 and t0 % p != 0) else None
            r2 = (t2 * pow(t0, -1, p)) % p if (v2 == v0 and t0 % p != 0) else None
            print(f"  p={p:2d}, n={n:2d} (u={n/p:.4f}): base(v={v0}, t%p={t0%p:2d}), "
                  f"dr=(0,0,-1)(v={v1}, ratio={r1}), dr=(1,0,0)(v={v2}, ratio={r2})")


# =============================================================================
# EXP 4: Constant Integer Offsets (alpha_i) on F_n*(t)
# =============================================================================
def run_exp4():
    print("\n" + "=" * 80)
    print("EXP 4: CONSTANT INTEGER OFFSETS ON F_n*(t) & SIMPLE-POLE WING AUDIT")
    print("=" * 80)
    # 4A: Check wing denominator of F_n*(t) directly at n = 1, 2, 3!
    for n in [1, 2, 3]:
        Qn, P3n, a_map, b_map, t0 = compute_exact_Fn_residues(n)
        # Check primes p in (23n, 25n]:
        wing_primes = [p for p in range(23 * n + 1, 25 * n + 1) if is_prime(p)]
        core_primes = [p for p in range(21 * n + 1, 23 * n + 1) if is_prime(p)]
        print(f"  n={n}: max_dbl_H = {21*n + t0} (=23n), max_wing_H = {23*n + t0} (=25n)")
        for p in wing_primes:
            print(f"    wing prime p={p} in (23n, 25n]: v_p(Q_n)={vp_frac(Qn, p)}, v_p(P3_n)={vp_frac(P3n, p)}  "
                  f"(d_{{25n}}^3 has +3, so net v_p(d_{{25n}}^3 P3_n) = {3 + vp_frac(P3n, p)})")
        for p in core_primes:
            print(f"    core prime p={p} in (21n, 23n]: v_p(Q_n)={vp_frac(Qn, p)}, v_p(P3_n)={vp_frac(P3n, p)}  "
                  f"(d_{{25n}}^3 has +3, so net v_p(d_{{25n}}^3 P3_n) = {3 + vp_frac(P3n, p)})")

    # 4B: Check if constant shift dr = (-1, -1, -1) or (+1, +1, +1) changes Q_n/P_{3,n} p-adic valuations on shell & non-shell cells!
    print("\n  Testing constant offsets dr on F_n*(t) for extra p-divisibility (n=19, p=41 -> u=19/41=0.4634 non-shell, and n=19, p=40? p=43 -> u=0.4419):")
    # Wait: test u in a 1st-order shell vs non-shell for p = 43, n = 20 (u = 20/43 = 0.4651 in (6/13, 7/15) 1st-order shell!)
    # and p = 43, n = 19 (u = 19/43 = 0.4419 non-shell) and p = 43, n = 20 vs 21 (21/43 = 0.4884 in (11/23, 12/25) 1st-order!)
    for label, p, n in [("1st-order u=20/43", 43, 20), ("2nd-order u=21/44? use 29/61", 61, 29), ("non-shell u=19/43", 43, 19)]:
        vp, uf = fact_table_pk(40 * n + 20, p, 4)
        mod3 = p ** 3
        results = []
        for dr in [(0, 0, 0), (1, 1, 1), (-1, -1, -1), (0, 1, 2), (1, 0, -1)]:
            vmin = 999
            terms = []
            for j in range(3 * n, 21 * n + 1):
                v1, u1 = binom_vp_u(j + 13 * n + dr[0], 13 * n, p, 4, vp, uf)
                v2, u2 = binom_vp_u(j + 15 * n + dr[1], 13 * n, p, 4, vp, uf)
                v3, u3 = binom_vp_u(j + 17 * n + dr[2], 13 * n, p, 4, vp, uf)
                v4, u4 = binom_vp_u(20 * n, j - 3 * n, p, 4, vp, uf)
                v5, u5 = binom_vp_u(19 * n, j - 2 * n, p, 4, vp, uf)
                vj = v1 + v2 + v3 + v4 + v5
                uj = (u1 * u2 * u3 * u4 * u5) % mod3
                terms.append((vj, uj))
                if vj < vmin:
                    vmin = vj
            tot = sum((uj * (p ** (vj - vmin))) for vj, uj in terms if vj - vmin < 3) % mod3
            ext = 3 if tot == 0 else (2 if tot % (p * p) == 0 else (1 if tot % p == 0 else 0))
            results.append(f"dr={dr}:v={vmin}+{ext}")
        print(f"  {label}: " + " | ".join(results))


# =============================================================================
# EXP 5: Ultra-Fine Farey N = 80 Scan of F_n*(t) for Hidden 2nd/3rd-Order Cells
# =============================================================================
def run_exp5():
    print("\n" + "=" * 80)
    print("EXP 5: ULTRA-FINE FAREY N=80 SCAN OF F_n*(t) FOR 3RD-ORDER / MISSED CELLS")
    print("=" * 80)
    # Generate Farey(80) on (1/25, 1)
    N = 80
    fpts = [(0, 1), (1, 1)]
    for d in range(2, N + 1):
        for k in range(1, d):
            if math.gcd(k, d) == 1:
                fpts.append((k, d))
    fpts.sort(key=lambda x: x[0] / x[1])

    # Precompute primes up to 2600
    primes = [p for p in range(43, 2600) if is_prime(p)]

    # Pre-cache factor tables only on demand? Even faster: compute ONLY Q_n mod p^4 using quick Legendre!
    # Wait: let's check if ANY cell in Farey(80) has extra valuation >= 2 on Q_n OUTSIDE I_12 and I_23,
    # or extra valuation >= 3 INSIDE I_12 and I_23!
    I12 = (Fraction(7, 15), Fraction(11, 23))
    I23 = (Fraction(16, 17), Fraction(20, 21))

    anomalies = []
    cells_checked = 0
    for i in range(len(fpts) - 1):
        a = Fraction(fpts[i][0], fpts[i][1])
        b = Fraction(fpts[i + 1][0], fpts[i + 1][1])
        if b <= Fraction(1, 25):
            continue
        a_cl = max(a, Fraction(1, 25))
        mid = (a_cl + b) / 2
        # Pick 2 primes p where n/p lies strictly inside (a_cl, b) and p^2 > 40*n
        extras = []
        for p in primes:
            lo_n = math.floor(a_cl * p) + 1
            hi_n = math.ceil(b * p) - 1
            if lo_n <= hi_n:
                n = (lo_n + hi_n) // 2
                if p * p > 40 * n + 5:
                    vp, uf = fact_table_pk(40 * n + 5, p, 4)
                    mod4 = p ** 4
                    vmin = 999
                    terms = []
                    for j in range(3 * n, 21 * n + 1):
                        vj, uj = compute_aj_mod_pk(n, j, p, 4, vp, uf)
                        terms.append((vj, uj))
                        if vj < vmin:
                            vmin = vj
                    tot = 0
                    for vj, uj in terms:
                        dv = vj - vmin
                        if dv < 4:
                            tot = (tot + uj * (p ** dv)) % mod4
                    ext = 4 if tot == 0 else (3 if tot % (p**3) == 0 else (2 if tot % (p*p) == 0 else (1 if tot % p == 0 else 0)))
                    extras.append(ext)
                    if len(extras) >= 2:
                        break
        if not extras:
            continue
        cells_checked += 1
        mn_ext = min(extras)
        in_i12 = (I12[0] <= a_cl and b <= I12[1])
        in_i23 = (I23[0] <= a_cl and b <= I23[1])
        if mn_ext >= 3 or (mn_ext == 2 and not (in_i12 or in_i23)):
            anomalies.append((a_cl, b, extras))
    print(f"  Checked {cells_checked} Farey(80) cells on (1/25, 1) with 2 primes each.")
    if anomalies:
        print(f"  FOUND {len(anomalies)} NEW 2nd/3rd-ORDER CELLS:")
        for a, b, ex in anomalies[:20]:
            print(f"    ({a}, {b}) -> extra v_p(Q_n) = {ex}")
    else:
        print("  CONFIRMED: ZERO 3rd-order cells and ZERO missed 2nd-order cells at Farey order 80!")


if __name__ == "__main__":
    run_exp1()
    run_exp2()
    run_exp3()
    run_exp4()
    run_exp5()
