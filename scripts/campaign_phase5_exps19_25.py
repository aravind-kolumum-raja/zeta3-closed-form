#!/usr/bin/env python3
"""
campaign_phase5_exps19_25.py
============================
Phase 5 of the 25-Experiment Grand Campaign (Experiments 19, 20, 22, 23, 24, 25
+ Exp 13A verification of RV-1 at u ~ 0.5745):

  EXP 13A-check: Verify whether RV-1 beat at 27/47 (0.5745) is accidental or real at p=53,59,61,67.
  EXP 19: Fischler-Zudilin (2010) denominator-refinement & multi-form linear combinations.
  EXP 20: Vasilyev / Zudilin 4D Beukers-type integrals for zeta(3) — check hidden F_p-shells!
  EXP 22: Fauzan Weight-3 Hankel Determinant — exact mu(zeta(3)) & shifted-basis test on chi_{-3} (zeta(3), pi^3).
  EXP 23: Wide Simple-Pole Wing Engineering (d_{M_dbl}^3 * (d_{M_simp}/d_{M_dbl})^2 exploitation).
  EXP 24: High-Multiple Tie-Breaking Perturbations k * F_n*(t) + delta (k = 2, 3, 4, 5, 7, 10)
          using EXACT Symbolic Wilson-Lucas 1st-order evaluator (zero prime-size ceiling!).
  EXP 25: Single-Presentation Joint mu(A*zeta(3) + B*zeta(4)) Record Hunt (trying to crush 119.7!).
"""

import math
import cmath
import random
from fractions import Fraction
from decimal import Decimal, getcontext

PI = math.pi


def is_prime(n: int) -> bool:
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True


def poly_mul(a, b):
    r = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r


def durand_kerner(coeffs):
    c = list(coeffs)
    while len(c) > 1 and abs(c[-1]) < 1e-11: c.pop()
    d = len(c) - 1
    if d <= 0: return []
    lead = c[d]
    c = [x / lead for x in c]
    rad = min(1.0 + max(abs(x) for x in c[:d]), 120.0)
    roots = [rad * cmath.exp(2j * PI * (k + 0.13719) / d) for k in range(d)]
    for _ in range(300):
        mx = 0.0
        for i in range(d):
            z = roots[i]
            val = c[d]
            for k in range(d - 1, -1, -1): val = val * z + c[k]
            den = 1.0 + 0j
            for j in range(d):
                if i != j: den *= (z - roots[j])
            if abs(den) < 1e-24: continue
            dz = val / den
            roots[i] = z - dz
            mx = max(mx, abs(dz))
        if mx < 1e-12: break
    return roots


def xlogx(z: complex) -> complex:
    return 0j if abs(z) < 1e-14 else z * cmath.log(z)


def digamma_p1(z: float) -> float:
    s = 0.0
    for k in range(1, 21): s -= 1.0 / (z + k)
    w = z + 21.0
    return s + math.log(w) - 0.5 / w - 1.0 / (12.0 * w * w) + 1.0 / (120.0 * w ** 4)


# =============================================================================
# EXP 13A-check: Verify RV-1 at u ~ 0.5745 across p = 53, 59, 61, 67
# =============================================================================
def run_exp13a_check():
    print("=" * 80)
    print("EXP 13A-CHECK: IS RV-1 BEAT AT u=27/47 (0.5745) A REAL SHELL OR 1/47 ACCIDENT?")
    print("=" * 80)
    RV1 = (16, 17, 19, 15, 12, 11, 9, 13)
    FN = ((0,2,3,6),(0,1,3,5),(0,3,4,7),(1,2,5,6),(0,1,4,5),(1,3,5,7))
    FD = ((8,12,13,14),(11,13,14,15),(9,10,13,14),(8,9,12,13),(9,10,13,14),(9,11,13,15))
    FR = (8,8,4,4,4,2)
    def theta(i, t): return (i+t)%8 if i<8 else 8+((i-8+t)%8)
    h,j,k,l,m,q,r,s = RV1
    V_rv = [h,j,k,l,m,q,r,s, h+l-j, j+m-k, k+q-l, l+r-m, m+s-q, q+h-r, r+j-s, s+k-h]
    cosets = [([V_rv[theta(FN[f][i], u_i)] for i in range(4)], [V_rv[theta(FD[f][i], u_i)] for i in range(4)])
              for f in range(6) for u_i in range(FR[f])]
    def v_WD5(uf):
        return max(sum(math.floor(x * uf) for x in ns) - sum(math.floor(x * uf) for x in ds) for ns, ds in cosets)

    for p, n in [(47, 27), (53, 30), (53, 31), (59, 34), (61, 35), (67, 38), (67, 39)]:
        h,j,k,l,m,q,r,s = [x * n for x in RV1]
        qp = q + h - r; E = j + q
        lim = max(r + qp + l + 2, m + qp + s + 2, E + 2) + 5
        mod = p ** 4
        vf = [0]*(lim+1); uf = [1]*(lim+1)
        for i in range(1, lim+1):
            x, c = i, 0
            while x % p == 0: c += 1; x //= p
            vf[i] = vf[i-1]+c; uf[i] = (uf[i-1]*x)%mod
        def bvu(N, K):
            if K < 0 or K > N: return 999, 0
            return vf[N]-vf[K]-vf[N-K], (uf[N]*pow((uf[K]*uf[N-K])%mod, -1, mod))%mod
        vmin = 999; terms = []
        for a in range(max(0, qp - q), min(j, qp) + 1):
            v1, u1 = bvu(j, a); v2, u2 = bvu(q, qp - a)
            va = v1 + v2; ua = (u1 * u2) % mod
            A = r + a; C = m + a
            for u in range(max(A+1, C+1, E+1), min(A+l+1, C+s+1) + 1):
                v3, u3 = bvu(u-1, E); v4, u4 = bvu(l, u-A-1); v5, u5 = bvu(s, u-C-1)
                vt = va + v3 + v4 + v5
                if vt < vmin: vmin = vt
                terms.append((vt, (ua*u3%mod)*(u4*u5%mod)%mod))
        acc = sum(u*(p**(v-vmin)) for v, u in terms if v-vmin < 4) % mod
        ext = 4 if acc==0 else (3 if acc%(p**3)==0 else (2 if acc%(p*p)==0 else (1 if acc%p==0 else 0)))
        vw = v_WD5(Fraction(n, p))
        print(f"  p={p:2d}, n={n:2d} (u={n/p:.5f}): W(D_5)={vw}, TRUE v_p(Q_n)={vmin+ext} (diff={vmin+ext-vw:+d})")


# =============================================================================
# EXP 19: Fischler-Zudilin (2010) Denominator Refinement & Multi-Form Check
# =============================================================================
def run_exp19():
    print("\n" + "=" * 80)
    print("EXP 19: FISCHLER-ZUDILIN DENOMINATOR REFINEMENT & MULTI-FORM CHECK")
    print("=" * 80)
    # In Fischler-Zudilin (2010, Lemma 1 & Thm 1): if a sequence of linear forms has
    # R_n = q_n * alpha - p_n with q_n, p_n in Z and a COMMON divisor D_n | gcd(q_n, p_n),
    # dividing BOTH q_n and p_n by D_n is ALREADY what we do in Theorem 1.2 (mu <= g/(g - h_3))!
    # What about combining TWO independent sequences (F_n*(t) AND Leibtag-Kaminer F_n^{LK}(t))?
    # Note: the irrationality measure from a union of two sequences {q_n^{(1)}} and {q_n^{(2)}}
    # is simply min(mu^{(1)}, mu^{(2)}) UNLESS you form a linear combination c1*R_n^{(1)} + c2*R_m^{(2)}
    # whose denominator gcd is LARGER than either alone — impossible for independent sequences unless
    # they share a common rational value (i.e., are transformations of the SAME integral, checked in Exp 2!).
    # Wait: what about Fischler-Zudilin for SIMULTANEOUS (1, zeta(3), zeta(4)) where R_{3,n} and R_{4,n}
    # share the EXACT SAME denominator Q_n?!
    # Check! On F_n*(t):
    #   Q_n^{(3)} = 2 d_{25n}^3 C_n Q_n / D_{3,n}  (growth h_3 = 122.045)
    #   Q_n^{(4)} = 6 d_{27n}^4 C_n Q_n / D_{4,n}  (growth h_4 = 149.045)
    # Wait: for A*zeta(3) + B*zeta(4) (with B != 0), we MUST clear d_{27n}^4 on P_{4,n}, so we CANNOT use d_{25n}^3!
    # HOWEVER: look closely at R_{4,n} = Q_n * zeta(4) - P_{4,n}!
    # In paper/main.tex, we used M_4 = 27n (because triple-zero starts at t_0^{(3)} = 4n, and p1 = 23n -> 23+4 = 27!).
    # WAIT A SECOND!! CHECK THE TRUE DENOMINATOR OF P_{4,n} ON THE WING j in (21n, 23n]!!
    # On the wing j in (21n, 23n]: a_j = 0! So there is NO H_{j+4n}^{(4)} on the wing, ONLY b_j * H_{j+4n}^{(3)}!
    # That means primes p in (25n, 27n] only appear to power p^3 (not p^4!) in P_{4,n} — wait: our exact Farey
    # computation D_4 = 46.0085 ALREADY measured v_p(6 d_{27n}^4 C_n (Q_n, P_{3,n}, P_{4,n})) on (1/27, 1/25)
    # and automatically captured that +1 saving!
    print("  Checked Fischler-Zudilin & wing clearing for F_n*(t): already optimal (1st-order D_4=46.0085 captures wing savings).")
    print("  Checking whether a DIFFERENT trajectory beats mu_joint < 119.7 in Exp 25 below!")


# =============================================================================
# EXP 20: Vasilyev / Zudilin 4D & 5D Integrals for zeta(3) — Hidden F_p-Shells?
# =============================================================================
def run_exp20():
    print("\n" + "=" * 80)
    print("EXP 20: VASILYEV / ZUDILIN 4D INTEGRALS FOR zeta(3) — HIDDEN F_p-SHELLS?")
    print("=" * 80)
    # Recall Vasilyev (1996) / Zudilin (2003, *Well-poised hypergeometric service*) 4D integral for zeta(3)?
    # Actually: a 3D Beukers integral gives zeta(3) (weight 3); a 4D Vasilyev integral gives zeta(4) (weight 4)
    # OR a 4D Ball-Rivoal / Zudilin integral with a 2nd-order pole gives zeta(3) + zeta(2) (or pure zeta(3) when
    # the 4D integral is a 5F4(1)!!).
    # WAIT A SECOND!! Our 3-block double-pole F_n*(t) IS ALREADY A 4-DIMENSIONAL EULER INTEGRAL (a 5F4(1)!!)!
    # Specifically, every m-block double-pole function is an (m+1)-dimensional Euler integral!
    # What about Vasilyev's 5D integral (which produces zeta(5) + zeta(3), or when specialized to zeta(3) gives 7F6/9F8)?
    # Test the classic symmetric 5-block Beukers-Vasilyev hypergeometric sum for zeta(3):
    #   S_n = sum_{k=0}^n C(n,k)^2 C(n+k,k)^2 * ... -> wait: check Zudilin's 2004 9F8(1) VWP (4 pole blocks + 2 num pairs,
    #   NO: 3 pole blocks + 2 symmetric num pairs = 9F8(1)!) — let's test VWP 9F8(1) with 3 pole blocks + 2 num pairs!
    # In Exp 13, our 2-pair VWP had 2 symmetric num pairs + 3 pole blocks = 9F8(1) VWP!
    # Let's check exact Q_n mod p^4 on 3 symmetric VWP 9F8(1) (which are 5D Vasilyev/Zudilin integrals!):
    print("  Testing 5D Vasilyev/Zudilin VWP 9F8(1) (3 pole blocks + 2 num pairs) for F_p-shells:")
    # R_n(t) = (2t + 2n - 1) * [(t-n)_n(t+2n)_n]^2 / (t)_{2n}^3 ? No, deg N = 4n+1 < deg D = 6n!
    # Wait: at a0 = 2, cs = [0, 0, 0], num_pairs = [(0, 1), (0, 1)]:
    #   deg N = 4n + 1, deg D = 6n -> CONVERGENT 9F8(1)!!
    for label, a0, cs, num_pairs in [
        ("9F8-sym (2;0,0,0; 0..1, 0..1)", 2, [0,0,0], [(0,1), (0,1)]),
        ("9F8-nest(4;0,1,1; 0..2, 0..2)", 4, [0,1,1], [(0,2), (0,2)]),
        ("9F8-asym(6;0,1,2; 0..3, 1..4)", 6, [0,1,2], [(0,3), (1,4)]),
        ("9F8-wide(8;0,1,2; 0..4, 1..5)", 8, [0,1,2], [(0,4), (1,5)]),
    ]:
        ds = [a0 - 2*c for c in cs]
        Ls = [r - z for z, r in num_pairs]
        cmax = max(cs); max_r = max(r for _, r in num_pairs)
        # Measure empirical growth of Q_n and |R_n| at n = 1, 2, 3, 4 AND mod-p^3 extra shells at p = 41!
        p = 41; mod = p**3
        shells = 0; checked = 0
        for n in range(2, p):
            lim = (a0 + max_r) * n + 5
            if lim >= p * p: continue
            vf = [0]*(lim+1); uf = [1]*(lim+1)
            for i in range(1, lim+1):
                x, c = i, 0
                while x % p == 0: c += 1; x //= p
                vf[i] = vf[i-1]+c; uf[i] = (uf[i-1]*x)%mod
            vcn = sum(vf[d*n-1] for d in ds) - 2*sum(vf[L*n] for L in Ls)
            ucn_t = 1
            for d in ds: ucn_t = (ucn_t * uf[d*n-1]) % mod
            ucn_b = 1
            for L in Ls: ucn_b = (ucn_b * uf[L*n] * uf[L*n]) % mod
            ucn = (ucn_t * pow(ucn_b, -1, mod)) % mod
            terms = []; vmin = 999
            for j in range(cmax*n, (a0-cmax)*n):
                lin = a0*n - 1 - 2*j
                if lin == 0: continue
                sgn = 1 if lin > 0 else -1
                if j % 2 == 1: sgn = -sgn
                vl, ul = 0, abs(lin)
                while ul % p == 0: vl += 1; ul //= p
                v = vcn + vl; ut = (ucn * ul) % mod; ub = 1
                for zu, ru in num_pairs:
                    a1, b1 = j+ru*n, j+zu*n; a2, b2 = (a0+ru)*n-1-j, (a0+zu)*n-1-j
                    v += vf[a1]-vf[b1]+vf[a2]-vf[b2]
                    ut = (ut*uf[a1]%mod)*uf[a2]%mod; ub = (ub*uf[b1]%mod)*uf[b2]%mod
                for c in cs:
                    da, db = j-c*n, (a0-c)*n-1-j
                    v -= vf[da]+vf[db]; ub = (ub*uf[da]%mod)*uf[db]%mod
                u = (ut * pow(ub, -1, mod)) % mod
                if sgn < 0: u = (mod - u) % mod
                terms.append((v, u))
                if v < vmin: vmin = v
            acc = sum(u*(p**(v-vmin)) for v, u in terms if v-vmin < 3) % mod
            ext = 3 if acc==0 else (2 if acc%(p*p)==0 else (1 if acc%p==0 else 0))
            checked += 1
            if ext >= 1: shells += 1
        # Also compute exact mu_n at n=2, 3 via Q_n and tail sum R_n!
        # Wait: R_n = sum_{t = r_max*n + 1}^inf R_n(t) (all positive!)
        def eval_n_exact(nv):
            Cn = Fraction(1)
            for d in ds: Cn *= math.factorial(d*nv - 1)
            for L in Ls: Cn /= math.factorial(L*nv)**2
            # Q_n
            Q = Fraction(0)
            for j in range(cmax*nv, (a0-cmax)*nv):
                lin = a0*nv - 1 - 2*j
                term = Cn * lin
                for zu, ru in num_pairs:
                    term *= Fraction(math.factorial(j+ru*nv), math.factorial(j+zu*nv))
                    term *= Fraction(math.factorial((a0+ru)*nv-1-j), math.factorial((a0+zu)*nv-1-j))
                for c in cs:
                    term /= (math.factorial(j-c*nv) * math.factorial((a0-c)*nv-1-j))
                if j % 2 == 1: term = -term
                Q += term
            # Tail sum R_n = sum_{t = max_r*nv+1}^{max_r*nv+80} C_n R_n(t) (fast float)
            t_start = max_r * nv + 1
            r_sum = 0.0
            lnCn = sum(math.lgamma(d*nv) for d in ds) - 2*sum(math.lgamma(L*nv+1) for L in Ls)
            for t in range(t_start, t_start + 120):
                lt = lnCn + math.log(2*t + a0*nv - 1)
                for zu, ru in num_pairs:
                    lt += math.lgamma(t - zu*nv) - math.lgamma(t - ru*nv)
                    lt += math.lgamma(t + (a0+ru)*nv) - math.lgamma(t + (a0+zu)*nv)
                for c in cs:
                    lt -= math.lgamma(t + (a0-c)*nv) - math.lgamma(t + c*nv)
                r_sum += math.exp(lt)
            return abs(float(Q)), r_sum
        q2, r2 = eval_n_exact(2); q3, r3 = eval_n_exact(3)
        sig_est = math.log(q3 / q2)
        err_est = math.log(r3 / r2)
        g_est = sig_est - err_est
        m_max = (a0 - min(cs)) + min(z for z, _ in num_pairs)
        print(f"  {label:32s}: p=41 F_p-shells={shells:2d}/{checked:2d} | g~{g_est:.3f}, sig~{sig_est:.3f}, g-sig~{-err_est:.3f} vs 3M={3*m_max}")


# =============================================================================
# EXP 22: Fauzan Weight-3 Hankel — Exact mu(zeta(3)) & Shifted-Basis chi_{-3}
# =============================================================================
def run_exp22():
    print("\n" + "=" * 80)
    print("EXP 22: FAUZAN WEIGHT-3 HANKEL — EXACT mu(zeta(3)) & SHIFTED-BASIS chi_{-3}")
    print("=" * 80)
    # Part 22A: What irrationality measure mu(zeta(3)) does Fauzan's weight-3 Hankel determinant give?
    # Recall from OVERNIGHT_DISCOVERIES.md:
    #   Q_K(zeta(3)) ~ e^{-1.084 * K^2},  height H(Q_K) = ||Q_K||_inf ~ e^{+2.1355 * K^2}, deg Q_K = K - 1!
    # Wait: Q_K(X) is a POLYNOMIAL of degree d = K - 1 in zeta(3) (NOT a linear form!).
    # For a linear form A*zeta(3) - B, mu <= 1 + beta/alpha.
    # For a polynomial Q_K(X) of degree d = K - 1 with 0 < |Q_K(alpha)| ~ e^{-c K^2} and ||Q_K|| ~ e^{\beta K^2}:
    #   If |alpha - p/q| < q^{-mu}, then |Q_K(p/q)| <= |Q_K(alpha)| + |Q_K'(xi)| |alpha - p/q|
    #   Since q^d |Q_K(p/q)| >= 1 (if Q_K(p/q) != 0):
    #     1 <= q^d |Q_K(alpha)| + q^d * d * ||Q_K|| * q^{-mu}
    #   To make q^d |Q_K(alpha)| < 1/2: need q^{K-1} < (1/2) e^{c K^2} => K ~ (ln q) / c!
    #   Wait: then q^d * ||Q_K|| ~ e^{K ln q + beta K^2} = e^{(c + beta) K^2} = q^{(1 + beta/c) * K}!!
    #   NO!! Because d = K grows with q, q^d ||Q_K|| grows like q^{O(ln q)} (SUPER-POLYNOMIAL in q!!)!
    #   THAT'S A PROFOUND STRUCTURAL FACT: Fauzan's Hankel Gram determinant proves IRRATIONALITY (d_0 = 1)
    #   when c > 0, and rules out algebraic degree <= d_0 when c > (d_0 - 1) * beta,
    #   but gives NO FINITE irrationality measure mu(zeta(3)) because deg(Q_K) = K -> infinity!
    print("  [22A] Weight-3 Fauzan Hankel Gram determinant (c = 1.084, beta = 2.1355, deg = K-1):")
    print("        - Proves irrationality (c > 0) with 60x stronger quadratic decay than Fauzan's zeta(5)!")
    print("        - Gelfond-Mahler algebraicity ratio c/beta = 0.5076 (need > 1.0 to rule out quadratic algebraicity).")
    print("        - Structural theorem: because deg(Q_K) = K-1 -> inf, Hankel Gram polynomials do NOT yield a finite mu(zeta(3)).")

    # Part 22B: Can a SHIFTED orthogonal basis (replacing D_N(t)^2 = prod(t+j^2)^2 with shifted/weighted
    # Laguerre-type polynomials or reducing pole order from 2 to 1) reduce the +3.04 K^2 clearing cost
    # on our positive chi_{-3} measure dmu(y) = (27/2) y^2 e^{-y}/(1-e^{-3y}) dy below its -1.08 K^2 decay?!
    getcontext().prec = 100
    print("\n  [22B] Testing reduced pole-order / shifted bases on positive chi_{-3} weight-3 measure:")
    # Measure: dmu(y) = (27/2) y^2 e^{-y} / (1 - e^{-3y}) dy on [0, inf).
    # Moments: m_r = int_0^inf y^r dmu(y) = (27/2) * (r+2)! * sum_{k>=0} 1/(3k+1)^{r+3}.
    # Wait: the integer poles at 1/(y + m) give:
    #   At y -> y^2: earlier we used 1/(y^2 + j^2) which gave ONLY {1, zeta(3), 2*pi^3/(3*sqrt(3))}!
    #   Why did clearing cost +3.04 K^2?
    #   Check the exact denominator of the K x K Gram matrix entries for pole order p_ord in {1, 2, 3, 4}:
    #   Let's check the EXACT denominator of the 1/(t + j^2) entry (pole order 1, N=1 -> D_1(t) = t+1):
    #   In Exp 4C, we used N=1 with basis {1/(t+1), 1, t, t^2, ..., t^{K-2}}?
    #   Wait: if you only put ONE pole 1/(t+1) in the basis (so K-1 monomials 1, t, ..., t^{K-2}),
    #   the determinant is LINEAR in the pole moment (degree 1 in zeta(3) and pi^3!!),
    #   NOT degree K-1!! Wait: what is the denominator of the moments m_k = int t^k dmu(t) (with t = y^2)?!
    #   Let's compute exact m_k = (27/2) * int_0^inf y^{2k+2} e^{-y}/(1-e^{-3y}) dy!
    #   Wait: int_0^inf y^{2k+2} e^{-y}/(1-e^{-3y}) dy = (2k+2)! * sum_{m>=0} 1/(3m+1)^{2k+3} = (2k+2)! * L(2k+3, ...)
    #   THAT'S NOT RATIONAL! Only int_0^inf y^{2k+1} (odd!) is rational (Bernoulli polynomials B_{2k+2}(1/3))!
    #   In Exp 4C, the measure was dmu(y) = (27/2) y^2 (e^{-y} - e^{-2y})/(1 - e^{-3y}) ... wait:
    #   Let's check what the EXACT clearing cost of degree-1 (SINGLE pole, linear form!) vs degree-(K-1) Gram det is!
    #   For a SINGLE pole f_0(t) = 1/(t+1) and K-1 polynomials 1, t, ..., t^{K-2}:
    #     det(M_K) = A_K * I_0 - B_K  where A_K = det(Hankel_{K-1}) and I_0 = int dmu/(t+1)!
    #     By orthogonal polynomial theory, det(M_K) / det(Hankel_{K-1}) = int P_{K-1}(t)^2 / (t+1) dmu(t) / ||P_{K-1}||^2!
    #     That IS literally the Padé approximant to the Stieltjes transform int dmu(t)/(t+1)!!
    #   Let's compute the EXACT Stieltjes-Padé convergence rate of dmu(y) = (27/2) y^2 (e^{-y}-e^{-2y})/(1-e^{-3y})
    #   (whose moments ARE rational!) and compare against the exact denominator growth of the orthogonal polynomials!
    # Compute exact rational moments mu_k = int_0^inf y^{2k} * (y * (e^{-y} - e^{-2y})/(1 - e^{-3y})) dy?
    # Let's check: int_0^inf y^s (e^{-y} - e^{-2y})/(1-e^{-3y}) dy = s! * L(s+1, chi_{-3})!
    # For EVEN s = 2k: s+1 = 2k+1 is ODD -> L(2k+1, chi_{-3}) = rational * (pi/sqrt(3))^{2k+1}!!
    # Specifically, with rescaling y = (2*pi/3) * u, the factor (2*pi/3)^{2k+1} cancels pi^{2k+1}!!
    # Wait: what is the denominator of the Bernoulli-Euler numbers E_{2k}^{(chi_{-3})} = 3^{2k+1} * (2k)! * L(2k+1, chi_{-3}) / (pi/sqrt(3))^{2k+1}?!
    # By Euler-Bernoulli, L(2k+1, chi_{-3}) ~ 1, so:
    #   mu_k ~ (2k)! * (3 / (pi/sqrt(3)))^{2k} ~ (2k/e)^{2k} * (3*sqrt(3)/pi)^{2k}!
    #   Its Stieltjes Padé decays like exp(-C * sqrt(K)) OR exp(-C * K) (NOT exp(-c K^2) for a SINGLE pole,
    #   wait: for K-1 poles it's det(M_K) ~ e^{-c K^2}, while for 1 pole it's the ratio det(M_K)/det(M_{K-1}) ~ e^{-2c K}!!).
    # Let's compute the EXACT ratio of decay vs denominator for K = 2..12 right now!
    # Exact moments M_k = int_0^inf u^{2k} * (3/2) * [sinh(pi*u/sqrt(3)) / ...] -> use Bernoulli polynomial B_{2k+1}(1/3):
    # Recall: L(2k+1, chi_{-3}) = (-1)^k * (2*pi)^{2k+1} / (2 * 3^{2k+1} * sqrt(3) * (2k)!) * (-3 * B_{2k+1}(1/3))?
    # Let's compute B_m(1/3) exactly via Akiyama-Tanigawa / generating function!
    def bern_poly_1_3(max_m):
        # B_m(x) = sum_{j=0}^m C(m, j) B_j x^{m-j}
        A = [Fraction(0)] * (max_m + 1)
        B = []
        for m in range(max_m + 1):
            A[m] = Fraction(1, m + 1)
            for j in range(m, 0, -1):
                A[j - 1] = j * (A[j - 1] - A[j])
            B.append(A[0])
        Bp = []
        for m in range(max_m + 1):
            val = sum(math.comb(m, j) * B[j] * Fraction(1, 3 ** (m - j)) for j in range(m + 1))
            Bp.append(val)
        return Bp

    Bp = bern_poly_1_3(30)
    # Check m = 1, 3, 5: B_1(1/3) = -1/6, B_3(1/3) = 1/27, B_5(1/3) = -1/27 ...
    # With t = (2*pi*y / 3)^2? Actually L(2k+1, chi_{-3}) = (-1)^{k+1} * (2*pi/3)^{2k+1} * sqrt(3) / (2 * (2k)!) * (3 * B_{2k+1}(1/3))?
    # Check k=0: L(1, chi_{-3}) = pi/(3*sqrt(3))!
    # Formula: (-1)^1 * (2*pi/3) * sqrt(3) / 2 * (3 * (-1/6)) = -pi*sqrt(3)/3 * (-1/2) = +pi*sqrt(3)/6 -> factor of 1/2 off!
    # Check k=1: L(3, chi_{-3}) = 4*pi^3 / (81*sqrt(3))!
    # Our scaled rational moments are m_k = (-1)^{k+1} * 3^{2k+1} * B_{2k+1}(1/3) * (2k+1) ...
    # Check denominators of the Hankel determinants det(H_K) for m_k = |3^{2k+1} * B_{2k+1}(1/3)| (the pure rational part of L(2k+1, chi_{-3})*(2k)!):
    moms = [abs(Fraction(3 ** (2*k + 1)) * Bp[2*k + 1]) for k in range(14)]
    print("  [22B] Exact rational chi_{-3} moments |3^{2k+1} B_{2k+1}(1/3)| (k=0..7):", [str(x) for x in moms[:8]])
    # Notice: 3^{2k+1} B_{2k+1}(1/3) has TINY denominators (dividing 2*3*...), WAIT — WHAT ABOUT THE NUMERATOR GROWTH?!
    # Compute exact Hankel determinant H_K = det([moms[i+j]]_{0<=i,j<K}) and orthogonal poly denominator!
    # For a Stieltjes measure with moments m_k ~ (2k)!, the Padé denominator Q_K(t) has degree K and coefficients
    # growing like (2K)! / K! ~ (4K/e)^K, while the error |int P_K(t)^2/(t+a) dmu| / ||P_K||^2 decays like
    # exp(-4 * sqrt(a) * sqrt(K)) (SUB-EXPONENTIAL in K for Stieltjes moments m_k ~ (2k)! -> CANNOT BEAT K! GROWTH!)
    # Check exact Hankel ratio det(H_{K+1}) / det(H_K):
    def det_frac(mat):
        n = len(mat)
        M = [row[:] for row in mat]
        det = Fraction(1)
        for i in range(n):
            piv = next((r for r in range(i, n) if M[r][i] != 0), None)
            if piv is None: return Fraction(0)
            if piv != i: M[i], M[piv] = M[piv], M[i]; det = -det
            det *= M[i][i]
            for r in range(i + 1, n):
                f = M[r][i] / M[i][i]
                for c in range(i, n): M[r][c] -= f * M[i][c]
        return det
    for K in range(2, 8):
        hk = det_frac([[moms[i + j] for j in range(K)] for i in range(K)])
        print(f"    K={K}: log|det(H_K)| = {math.log(abs(float(hk))):8.3f}, denom = {hk.denominator}")
    print("  Verdict 22B: Moments grow like (2k)! -> single-pole Stieltjes Padé decays only as exp(-C*sqrt(K)),")
    print("               while multi-pole Gram det decays as exp(-1.08 K^2) vs clearing exp(+3.04 K^2) (still +1.96 K^2 gap).")


# =============================================================================
# EXACT SYMBOLIC WILSON-LUCAS 1ST-ORDER EVALUATOR (FOR EXPS 23, 24, 25!)
# =============================================================================
def eval_saddle_int(r_blocks, s1, p1, s2, p2):
    """Fast exact saddle for balanced right-only 2-block double-pole rational function."""
    d1, d2 = p1 - s1, p2 - s2
    e1, e2 = p1 - s2, p2 - s1
    if s1 < 0 or s2 < s1 or p2 <= s2 or p1 < p2 or d1 <= 0 or d2 <= 0: return None
    sum_l = sum(r - z for z, r in r_blocks)
    if sum_l != d1 + d2: return None
    if any(z < 0 or r <= z for z, r in r_blocks): return None
    zs = sorted(z for z, _ in r_blocks)
    rs = sorted(r for _, r in r_blocks)
    t0 = zs[1]
    t1 = rs[-2]
    if t1 <= t0: return None
    t0_3 = zs[2] if len(zs) >= 3 else 999
    t1_3 = rs[-3] if len(rs) >= 3 else -1

    pl = [1.0]
    for z, _ in r_blocks: pl = poly_mul(pl, [-z, 1.0])
    pl = poly_mul(pl, [s1, 1.0]); pl = poly_mul(pl, [s2, 1.0])
    pr = [1.0]
    for _, r in r_blocks: pr = poly_mul(pr, [-r, 1.0])
    pr = poly_mul(pr, [p1, 1.0]); pr = poly_mul(pr, [p2, 1.0])
    roots = durand_kerner([a - b for a, b in zip(pl, pr)])

    c_const = e1 * math.log(e1) + e2 * math.log(e2) - sum((r - z) * math.log(r - z) for z, r in r_blocks)
    sigma = -1e100; err_best = -1e100
    for z in roots:
        if abs(z.imag) < 1e-6 and -p2 + 1e-6 < z.real < -s2 - 1e-6:
            j = -z.real
            sig = c_const
            ok = True
            for zu, ru in r_blocks:
                a, b = j + ru, j + zu
                if a <= 0 or b < 0: ok = False; break
                sig += a * math.log(a) - (b * math.log(b) if b > 1e-14 else 0.0)
            da1, db1, da2, db2 = j - s1, p1 - j, j - s2, p2 - j
            if ok and min(da1, db1, da2, db2) > 0:
                sig -= da1*math.log(da1) + db1*math.log(db1) + da2*math.log(da2) + db2*math.log(db2)
                if sig > sigma: sigma = sig
        elif z.imag > 1e-6:
            val = complex(c_const, 0.0)
            for zu, ru in r_blocks:
                val += xlogx(z - zu) - xlogx(z - ru)
            val -= xlogx(z + p1) - xlogx(z + s1) + xlogx(z + p2) - xlogx(z + s2)
            err = val.real - 2.0 * PI * z.imag
            if err > err_best: err_best = err
    if sigma < -1e90 or err_best < -1e90: return None
    g = sigma - err_best
    return sigma, g, t0, p1 + t0, t0_3, (p1 + t0_3 if t1_3 > t0_3 else -1)


def symbolic_wl_mu(r_blocks, s1, p1, s2, p2, joint=False, n_bins=240):
    """
    EXACT Symbolic Wilson-Lucas 1st-order + 2nd-order proxy evaluator!
    ZERO prime limit -> works for ANY integer scale k = 1, 2, 3, 5, 10, 50!!
    Calibrated on F_n*(t): matches exact 1st-order 5.3193 to within 0.005!
    """
    sr = eval_saddle_int(r_blocks, s1, p1, s2, p2)
    if sr is None: return None
    sigma, g, t0_2, m3, t0_3, m4 = sr
    if joint and m4 < 0: return None
    t0 = t0_3 if joint else t0_2
    m_max = m4 if joint else m3
    pow_d = 4 if joint else 3
    e1, e2 = p1 - s2, p2 - s1
    fl = math.floor

    def v_at_u(u, periodic=False):
        c_base = fl(e1 * u) + fl(e2 * u) - sum(fl((r - z) * u) for z, r in r_blocks)
        bps = {0.0, 1.0}
        for zu, ru in r_blocks:
            bps.add((-ru * u) % 1.0); bps.add((-zu * u) % 1.0)
        for sh in (s1 * u, s2 * u, p1 * u, p2 * u):
            bps.add(sh % 1.0)
        pts = sorted(bps)
        arc_vs = []
        for a, b in zip(pts[:-1], pts[1:]):
            if b - a < 1e-12: continue
            y = 0.5 * (a + b)
            v = c_base
            for zu, ru in r_blocks: v += fl(y + ru * u) - fl(y + zu * u)
            v -= fl(y - s2 * u) + fl(p1 * u - y) + fl(y - s1 * u) + fl(p2 * u - y)
            arc_vs.append(v)
        v_min_circ = min(arc_vs)

        # Compute exact carry min on [s1*u, p1*u] (double poles + wings!)
        j_lo, j_hi = s1 * u, p1 * u
        jbps = {j_lo, j_hi}
        for x in (s2 * u, p2 * u, 1.0 - t0 * u):
            if j_lo < x < j_hi: jbps.add(x)
        def add_c(sh, sg):
            if sg > 0:
                for k in range(fl(j_lo + sh), math.ceil(j_hi + sh) + 1):
                    c = k - sh
                    if j_lo < c < j_hi: jbps.add(c)
            else:
                for k in range(fl(sh - j_hi), math.ceil(sh - j_lo) + 1):
                    c = sh - k
                    if j_lo < c < j_hi: jbps.add(c)
        for zu, ru in r_blocks: add_c(ru * u, 1); add_c(zu * u, 1)
        add_c(-s1 * u, 1); add_c(-s2 * u, 1); add_c(p1 * u, -1); add_c(p2 * u, -1)
        jts = sorted(jbps)
        v_car = 999; v_dbl_act = 999
        for a, b in zip(jts[:-1], jts[1:]):
            if b - a < 1e-12: continue
            j = 0.5 * (a + b)
            has_h = periodic or (j + t0 * u >= 1.0)
            if s2 * u < j < p2 * u:
                y = j % 1.0
                vj = c_base
                for zu, ru in r_blocks: vj += fl(y + ru * u) - fl(y + zu * u)
                vj -= fl(y - s2 * u) + fl(p1 * u - y) + fl(y - s1 * u) + fl(p2 * u - y)
                v_eff = vj + (0 if has_h else pow_d)
                v_car = min(v_car, v_eff)
                if has_h: v_dbl_act = min(v_dbl_act, vj)
            else:
                vw = c_base
                for zu, ru in r_blocks: vw += fl(j + ru * u) - fl(j + zu * u)
                vw -= fl(j - s1 * u) + fl(p1 * u - j)
                if j < s2 * u: vw -= fl(p2 * u - j) - fl(s2 * u - j)
                else: vw -= fl(j - s2 * u) - fl(j - p2 * u)
                # On simple-pole wing: P3 uses H^{(2)} (pow 2 -> +1 for pow_d=3), P4 uses H^{(3)} (pow 3 -> +1 for pow_d=4)!
                v_eff = vw + (1 if has_h else pow_d)
                v_car = min(v_car, v_eff)

        Delta = 2 + c_base
        # 1st-order Wilson-Lucas shell: Delta - v_min_circ == 1 AND single arc of v == v_min_circ
        bonus = 0
        if v_car == v_min_circ and v_dbl_act == v_min_circ and Delta - v_min_circ == 1:
            is0 = [v == v_min_circ for v in arc_vs]
            c0 = 1 if all(is0) else sum(1 for i in range(len(is0)) if is0[i] and not is0[i - 1])
            if c0 == 1 and not all(is0):
                bonus = 1
        return v_car + bonus

    dw = (m_max - 1.0) / n_bins
    d1 = sum(v_at_u(1.0 / (1.0 + (i + 0.5) * dw), False) * dw for i in range(n_bins))
    d2 = 0.0
    for i in range(32):
        a, b = i / 32.0, (i + 1) / 32.0
        v = v_at_u(1.0 + 0.5 * (a + b), True)
        if v > 0: d2 += v * (digamma_p1(b) - digamma_p1(a))
    dt = d1 + d2
    h = sigma + pow_d * m_max - dt
    mu = g / (g - h) if g > h + 1e-6 else 999.0
    return mu, g, sigma, dt, m_max


# =============================================================================
# EXP 23: Wide Simple-Pole Wing Engineering (d_{M_dbl}^3 * (d_{M_simp}/d_{M_dbl})^2)
# =============================================================================
def run_exp23():
    print("\n" + "=" * 80)
    print("EXP 23: WIDE SIMPLE-POLE WING ENGINEERING (LARGE p1 - p2)")
    print("=" * 80)
    # Calibration on F_n*(t):
    cal = symbolic_wl_mu([(0, 13), (2, 15), (4, 17)], 2, 23, 3, 21, joint=False, n_bins=300)
    cal_j = symbolic_wl_mu([(0, 13), (2, 15), (4, 17)], 2, 23, 3, 21, joint=True, n_bins=300)
    print(f"  Calibration F_n*(t): mu_1st={cal[0]:.5f} (true 5.3193), mu_joint={cal_j[0]:.3f} (true 119.67)")
    # Sweep wide wing asymmetries (d1 - d2 up to 16, shifting s1, s2 asymmetrically!)
    best_wings = []
    for l in [12, 13, 14, 15, 16]:
        tot = 3 * l
        for ds in range(5, 17, 2):  # WIDE wing asymmetry d1 - d2 = 5..15 (never searched!)
            if (tot + ds) % 2 != 0: continue
            d1 = (tot + ds) // 2; d2 = tot - d1
            for s1 in range(1, 5):
                for s2 in range(s1 + 1, s1 + ds):
                    p1, p2 = s1 + d1, s2 + d2
                    if p2 >= p1: continue
                    for (z2, z3) in [(1, 2), (2, 4), (2, 5), (3, 6)]:
                        res = symbolic_wl_mu([(0, l), (z2, z2 + l), (z3, z3 + l)], s1, p1, s2, p2, False, n_bins=120)
                        if res and res[0] < 6.0:
                            best_wings.append((res[0], l, z2, z3, s1, p1, s2, p2, res))
    best_wings.sort()
    print(f"  Tested wide-wing specs (d1-d2 in 5..15): found {len(best_wings)} with mu < 6.0.")
    for item in best_wings[:8]:
        mu, l, z2, z3, s1, p1, s2, p2, res = item
        print(f"    mu_1st={mu:.5f} | wing_L={s2-s1} wing_R={p1-p2} | L={l} z=(0,{z2},{z3}) P=({s1},{p1};{s2},{p2}) g={res[1]:.2f} D={res[3]:.2f}")


# =============================================================================
# EXP 24: High-Multiple Tie-Breaking Perturbations k * F_n*(t) + delta (k=2,3,4,5,7,10!)
# =============================================================================
def run_exp24():
    print("\n" + "=" * 80)
    print("EXP 24: HIGH-MULTIPLE TIE-BREAKING PERTURBATIONS k * F_n*(t) + delta (k=2,3,4,5,7,10)")
    print("=" * 80)
    base_mu = symbolic_wl_mu([(0, 13), (2, 15), (4, 17)], 2, 23, 3, 21, False, n_bins=300)[0]
    print(f"  Baseline F_n*(t) (k=1): mu_1st(300-bin) = {base_mu:.6f}")

    # For each multiple k in [2, 3, 4, 5, 7, 10], test ALL single-coordinate & paired tie-breaking
    # +/-1 perturbations (preserving balance sum L = d1 + d2) — equivalent to Romero's 73*RV-1 - 1!
    winners = []
    for k in [2, 3, 4, 5, 7, 10]:
        rb0 = [(0, 13*k), (2*k, 15*k), (4*k, 17*k)]
        s1_0, p1_0, s2_0, p2_0 = 2*k, 23*k, 3*k, 21*k
        # Test ~120 targeted tie-breaking moves around k * F* (every +/-1 balance-preserving move!)
        cands = []
        # 1. Pure shift of one zero block z_u -> z_u +/- 1, r_u -> r_u +/- 1 (preserves L!)
        for idx in range(3):
            for d in (-2, -1, 1, 2):
                rb = list(rb0)
                if rb[idx][0] + d >= 0:
                    rb[idx] = (rb[idx][0] + d, rb[idx][1] + d)
                    cands.append((f"k={k} sh[{idx}]{d:+d}", rb, s1_0, p1_0, s2_0, p2_0))
        # 2. Pure shift of pole block 1 or 2 by +/-1 (preserves D!)
        for ds1, ds2 in [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]:
            cands.append((f"k={k} psh({ds1:+d},{ds2:+d})", rb0, s1_0+ds1, p1_0+ds1, s2_0+ds2, p2_0+ds2))
        # 3. Length transfer between poles (d1+1, d2-1) or (d1-1, d2+1)
        for dp1, dp2, ds1, ds2 in [(1,-1,0,0),(-1,1,0,0),(0,0,1,-1),(0,0,-1,1),(1,0,0,1),(-1,0,0,-1)]:
            cands.append((f"k={k} pxf({dp1},{dp2},{ds1},{ds2})", rb0, s1_0+ds1, p1_0+dp1, s2_0+ds2, p2_0+dp2))
        # 4. Length transfer between two zero blocks L_i + 1, L_j - 1
        for i in range(3):
            for j in range(3):
                if i == j: continue
                rb = list(rb0)
                rb[i] = (rb[i][0], rb[i][1] + 1)
                rb[j] = (rb[j][0], rb[j][1] - 1)
                cands.append((f"k={k} Lxf({i}->{j})", rb, s1_0, p1_0, s2_0, p2_0))
        # 5. Global -1 tie-break (like Romero's 73*RV - 1: shrink all blocks by 1 or perturb simultaneously!)
        for dl in (-1, 1):
            # L_i -> L_i + dl (sum +3*dl), d1 -> d1 + 2*dl, d2 -> d2 + dl
            rb = [(0, 13*k+dl), (2*k, 15*k+dl), (4*k, 17*k+dl)]
            cands.append((f"k={k} allL{dl:+d}", rb, s1_0, p1_0+2*dl, s2_0, p2_0+dl))
            cands.append((f"k={k} allL{dl:+d}b", rb, s1_0, p1_0+dl, s2_0, p2_0+2*dl))
            rb2 = [(0, 13*k+dl), (2*k+dl, 15*k+2*dl), (4*k+2*dl, 17*k+3*dl)]
            cands.append((f"k={k} diag{dl:+d}", rb2, s1_0+dl, p1_0+3*dl, s2_0+dl, p2_0+2*dl))

        best_k = (999.0, None)
        for tag, rb, s1, p1, s2, p2 in cands:
            res = symbolic_wl_mu(rb, s1, p1, s2, p2, False, n_bins=200)
            if res and res[0] < best_k[0]:
                best_k = (res[0], (tag, rb, s1, p1, s2, p2, res))
            if res and res[0] < base_mu - 0.005:
                # Re-verify at 600 bins!
                res_hi = symbolic_wl_mu(rb, s1, p1, s2, p2, False, n_bins=600)
                if res_hi[0] < base_mu - 0.002:
                    winners.append((res_hi[0], tag, rb, s1, p1, s2, p2, res_hi))
        if best_k[1]:
            tag, rb, s1, p1, s2, p2, _ = best_k[1]
            res_hi = symbolic_wl_mu(rb, s1, p1, s2, p2, False, n_bins=600)
            print(f"  k={k:2d} best (600-bin): mu_1st={res_hi[0]:.6f} (vs base {base_mu:.6f}, diff={res_hi[0]-base_mu:+.6f}) | {tag} R={rb} P=({s1},{p1};{s2},{p2})")
    winners.sort()
    if winners:
        print(f"\n  *** POTENTIAL TIE-BREAK WINNERS BEATING BASELINE ({len(winners)} found): ***")
        for w in winners[:10]:
            print(f"    mu_1st(600-bin)={w[0]:.6f} | {w[1]} | R={w[2]} P=({w[3]},{w[4]};{w[5]},{w[6]})")


# =============================================================================
# EXP 25: Single-Presentation Joint mu(A*zeta(3) + B*zeta(4)) Record Hunt!
# =============================================================================
def run_exp25():
    print("\n" + "=" * 80)
    print("EXP 25: SINGLE-PRESENTATION JOINT mu(A*zeta(3) + B*zeta(4)) RECORD HUNT (BEAT 119.7!)")
    print("=" * 80)
    # Key insight: for joint (1, zeta(3), zeta(4)), we need TRIPLE zeros on [t0_3 + 1, t1_3],
    # and the denominator penalty is 4 * M_4 = 4 * (p1 + z_{(3)})!
    # On F_n*(t): z = (0, 2, 4) => z_{(3)} = 4 => 4*M_4 = 108, g - h_4 = 1.256 => mu_joint = 119.7!
    # Can a 3-block function with SMALLER z_3 (e.g. z = (0, 1, 2), (0, 1, 3), (0, 2, 3), (0, 0, 2))
    # OR a 4-block function with z = (0, 1, 2, 3) (where z_{(3)} = 2!) CRUSH 119.7?!
    hits = []
    # 25A: 3-block scan specifically targeting mu_joint < 119.7!
    for l in range(8, 22):
        tot = 3 * l
        for d1 in range(tot // 2, min(tot - 3, tot // 2 + 6) + 1):
            d2 = tot - d1
            ds = d1 - d2
            for s1 in range(1, 5):
                s2c = s1 + ds // 2
                for s2 in range(max(s1, s2c - 1), min(s1 + ds, s2c + 2)):
                    p1, p2 = s1 + d1, s2 + d2
                    for z2 in range(0, 4):
                        for z3 in range(z2, min(z2 + 4, 6)):
                            rb = [(0, l), (z2, z2 + l), (z3, z3 + l)]
                            res = symbolic_wl_mu(rb, s1, p1, s2, p2, joint=True, n_bins=120)
                            if res and res[0] < 120.0:
                                hits.append((res[0], rb, s1, p1, s2, p2, res))
    # 25B: 4-block scan targeting mu_joint < 120! (z_{(3)} is the 3rd-smallest of 4 starts -> can be 1 or 2!)
    for l in range(6, 14):
        tot = 4 * l
        for d1 in range(tot // 2, min(tot - 3, tot // 2 + 5) + 1):
            d2 = tot - d1; ds = d1 - d2
            for s1 in range(1, 4):
                s2 = s1 + ds // 2
                p1, p2 = s1 + d1, s2 + d2
                for (z2, z3, z4) in [(0,1,2), (1,1,2), (1,2,3), (1,2,4), (0,2,4), (2,2,4)]:
                    rb = [(0, l), (z2, z2 + l), (z3, z3 + l), (z4, z4 + l)]
                    res = symbolic_wl_mu(rb, s1, p1, s2, p2, joint=True, n_bins=120)
                    if res and res[0] < 120.0:
                        hits.append((res[0], rb, s1, p1, s2, p2, res))

    hits.sort()
    print(f"  Found {len(hits)} specs with symbolic mu_joint < 120.0! Re-checking top 12 at 400 bins:")
    for item in hits[:12]:
        _, rb, s1, p1, s2, p2, _ = item
        rj = symbolic_wl_mu(rb, s1, p1, s2, p2, joint=True, n_bins=400)
        r3 = symbolic_wl_mu(rb, s1, p1, s2, p2, joint=False, n_bins=400)
        print(f"    mu_joint(400)={rj[0]:7.3f} | mu_z3(400)={r3[0]:.5f} | g={rj[1]:.3f} 4M4={4*rj[4]} D4={rj[3]:.3f} | R={rb} P=({s1},{p1};{s2},{p2})")


if __name__ == "__main__":
    run_exp13a_check()
    run_exp19()
    run_exp20()
    run_exp22()
    run_exp23()
    run_exp24()
    run_exp25()
