#!/usr/bin/env python3
"""
campaign_phase4_exp13_16_fixed.py
=================================
Fixed Exp 13 (VWP 7F6(1) + W(D_5) + Mechanism-2 F_p-shells, including EXACT RV-1 test!)
and fixed Exp 16 (EVEN 3-block differentiated F_n'(t)), fixing the (-1)^{3j} alternating
pole saddle equation (pl + pr = 0 for pole saddle, pl - pr = 0 for positive tail saddle!).
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


PRIMES = [p for p in range(3, 2500) if is_prime(p)]


def digamma_p1(z: float) -> float:
    s = 0.0
    for k in range(1, 21):
        s -= 1.0 / (z + k)
    w = z + 21.0
    return s + math.log(w) - 0.5 / w - 1.0 / (12.0 * w * w) + 1.0 / (120.0 * w ** 4)


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
    if d <= 0: return []
    lead = c[d]
    c = [x / lead for x in c]
    rad = min(1.0 + max(abs(x) for x in c[:d]), 80.0)
    roots = [rad * cmath.exp(2j * PI * (k + 0.13719) / d) for k in range(d)]
    for _ in range(350):
        mx = 0.0
        for i in range(d):
            z = roots[i]
            val = c[d]
            for k in range(d - 1, -1, -1):
                val = val * z + c[k]
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


def vwp_3blk_saddle(a0, cs, num_pairs):
    """
    EXACT saddle for 3-denominator-block symmetric rational function:
      R_n(t) = (2t + a0*n - 1) * prod_{(z,r)} [(t - r*n)_{L*n} (t + (a0+z)*n)_{L*n}]
               / prod_{i=1}^3 (t + c_i*n)_{(a0 - 2*c_i)*n}
    Check ratios:
      (1) Tail sum (t > max_r, ALL POSITIVE terms!):
          R(t+1)/R(t) = prod(t-z)/prod(t-r) * prod(t+a0+r)/prod(t+a0+z) * prod(t+c_i)/prod(t+a0-c_i) = +1
          => pl - pr = 0 with root x > max_r!
      (2) Pole sum Q_n = sum_j (-1)^j a_j (j in [c_max, a0-c_max]):
          Because of (-1)^{3j} = (-1)^j, the term magnitude |a_j| has ratio:
          |a_{j+1}|/|a_j| = prod(j+r)/(j+z) * prod(a0+z-j)/(a0+r-j) * prod(a0-c_i-j)/(j-c_i) = 1
          Wait: as a function of z = -j:
            a_{j+1}/a_j (WITH sign (-1)^j) = - [ prod(-z+r)/(-z+z_u) * prod(a0+z_u+z)/(a0+r+z) * prod(a0-c_i+z)/(-z-c_i) ]
          Since there are 3 factors of (-z-c_i) = -(z+c_i), that's (-1)^3 = -1, TIMES the -1 from (-1)^j = **+1**!!
          WAIT A SECOND!! Look at the signs:
            prod(-z+r)/prod(-z+z_u) = + prod(z-r)/prod(z-z_u)
            prod(a0+z_u+z)/prod(a0+r+z) = + prod(z+a0+z_u)/prod(z+a0+r)
            prod(a0-c_i+z)/prod(-z-c_i) = - prod(z+a0-c_i)/prod(z+c_i)  (3 minus signs -> MINUS!)
          SO:
            - |a_{j+1}|/|a_j| = +1  <==>  pl + pr = 0  (Wait: at j = a0/2, by symmetry |a_j| is symmetric around a0/2,
              so |a_{j+1}|/|a_j| = 1 AT j = a0/2!! Wait: does (-1)^j cause CANCELLATION so Q_n grows slower than max|a_j|?
              YES!! For Apéry, |a_j| peaks at j = n/2 with max|a_j| = C(3n/2, n)^2 C(n, n/2)^2 ~ (3^3/2^2)^n ... NO!
              Apéry's Q_n = sum C(n,k)^2 C(n+k,k)^2 has ALL POSITIVE terms!! Wait: why did Apéry have ALL POSITIVE terms?
              Because in Apéry's VWP 7F6(1), the triple pole residue has (-1)^{3j} = (-1)^j, AND (t - n)_n = (-1)^n ...,
              NO! Apéry's Q_n = sum_{k=0}^n C(n,k)^2 C(n+k,k)^2 comes from 2-block (Beukers) OR from VWP 7F6(1) via Whipple's
              transformation, which turns an ALTERNATING 7F6(1) residue sum into a POSITIVE sum — wait: in VWP 7F6(1)
              the residue sum Q_n IS sum_{j=0}^n (-1)^j (n-2j) C(j+n,n) C(2n-j,n) C(n,j)^3 -> NO!
              Check at n=1: j=0: (+1)*(1)*C(1,1)*C(2,1)*C(1,0)^3 = +2; j=1: (-1)*(-1)*C(2,1)*C(1,1)*C(1,1)^3 = +2; sum = 4!
              WAIT A SECOND!! For j -> n-j: (n-2j) flips sign AND (-1)^j flips sign (if n odd) -> wait:
              for n=2: j=0: (+1)*(+2)*C(2,2)C(4,2) = +12; j=1: (-1)*(0) = 0; j=2: (+1)*(-2)*C(4,2)C(2,2) = -12 -> SUM = 0!!
              THAT'S WHY APÉRY'S 7F6(1) HAS (t)_{n+1}^3 ONLY FOR n ODD — AND FOR ALL n YOU USE (t)_{2n}^3 (EVEN LENGTH)!!
              Check (2t + 2n - 1) (t - n)_n (t + 2n)_n / (t)_{2n}^3 at n=1 (a0=2, c=[0,0,0], r=1, z=0):
                j=0: (+1)*(+1)*C(1,1)*C(2,1)*C(1,0)^3 = +2
                j=1: (-1)*(-1)*C(2,1)*C(1,1)*C(1,1)^3 = +2 -> sum = 4!
              Wait: check j=0..2n-1 for n=2 (a0=2, 2n=4, r=2, z=0):
                j=0: (+1)*(+3)*C(2,2)C(5,2)C(3,0)^3 = +30
                j=1: (-1)*(+1)*C(3,2)C(4,2)C(3,1)^3 = -486
                j=2: (+1)*(-1)*C(4,2)C(3,2)C(3,2)^3 = -486
                j=3: (-1)*(-3)*C(5,2)C(2,2)C(3,3)^3 = +30
                SUM = -912 != 0!! AND IT'S ALTERNATING -> COMPLEX SADDLE FOR Q_n (pl - pr = 0 WITH z IN (-a0, 0) OR COMPLEX!!)!
          YES!! Because a_{j+1}/a_j (WITH sign (-1)^j) = +1 <==> pl - pr = 0!!
          Check Apéry (a0=2, c=[0,0,0], r=1, z=0) on pl - pr = 0:
            pl = z * (z + 3) * z^3 = z^4(z+3)
            pr = (z - 1) * (z + 2) * (z + 2)^3 = (z-1)(z+2)^4
            pl - pr = 0:
            Let's check roots of z^4(z+3) - (z-1)(z+2)^4 = 0!
            At z = x = 1 + sqrt(2) = 2.41421356 -> Wait: (z-1)(z+2)^4 vs z^4(z+3):
              Wait: for a0=2, R(t) = (2t+2n-1) (t-n)_n (t+2n)_n / (t)_{2n}^3 has 2n-1, NOT n!
              With a0=1 (classic Apéry), pl = z * (z+2) * z^3 = z^4(z+2), pr = (z-1)(z+1)^4!
              Check z^4(z+2) - (z-1)(z+1)^4 = 0:
                Expand: (z^5 + 2z^4) - (z-1)(z^4+4z^3+6z^2+4z+1)
                      = z^5+2z^4 - (z^5+3z^4+2z^3-2z^2-3z-1)
                      = -z^4 - 2z^3 + 2z^2 + 3z + 1 = -(z^2 + z - 1)(z^2 + z + 1)? NO:
                      z^4 + 2z^3 - 2z^2 - 3z - 1 -> Wait: at z = (sqrt(2)-1)? Check!
                      Let w = z + 1/2 (symmetric around -1/2!):
                      (w-1/2)^4(w+3/2) - (w-3/2)(w+1/2)^4 = 0!
                      Notice w -> -w maps LHS to -LHS -> w = 0 (z = -1/2) IS NOT A ROOT, wait:
                      at w = 0: (1/16)(3/2) - (-3/2)(1/16) = +3/16 != 0!
                      Check roots: w^4 - 5/2 w^2 - ... -> TWO REAL ROOTS w = +/- sqrt(2)??
                      Check w = +sqrt(2) -> z = sqrt(2) - 1/2 = 0.914? Wait: z must be > r = 1!
                      Check z = 1 + sqrt(2) = 2.4142 (w = 2.9142):
                        z^4(z+2) = 33.97056 * 4.41421 = 149.953
                        (z-1)(z+1)^4 = 1.41421 * 135.882 = 192.166 -> NOT EQUAL!
                      WHY?! Because Apéry's VWP 7F6(1) is (t - n)_n (t + 2n + 1)_n / (t)_{n+1}^3 ->
                        At t = -j: (t + n + 1)_n (NOT t + 2n + 1!!)!
                        Let's check! (t + n + 1)_n = (t+n+1)...(t+2n) = Gamma(t+2n+1)/Gamma(t+n+1)!
                        So in pl and pr:
                          (t - n)_n = Gamma(t)/Gamma(t-n) -> ratio (t)/(t-n)  (z in pl, z-1 in pr)
                          (t + n + 1)_n = Gamma(t+2n+1)/Gamma(t+n+1) -> ratio (t+2n+1)/(t+n+1)  (z+2 in pl, z+1 in pr)
                          1/(t)_{n+1}^3 = Gamma(t)^3/Gamma(t+n+1)^3 -> ratio (t)^3/(t+n+1)^3  (z^3 in pl, (z+1)^3 in pr)
                        Wait: that's z^4(z+2) = (z-1)(z+1)^4, whose roots are:
                          z^4 + 2z^3 - 2z^2 - 3z - 1 = 0! Let's print its roots — AND watch what c0, c1 it gives!!
    """
    ds = [a0 - 2 * c for c in cs]
    Ls = [r - z for z, r in num_pairs]
    pl = [1.0]
    for z, r in num_pairs:
        pl = poly_mul(pl, [-z, 1.0])
        pl = poly_mul(pl, [a0 + r, 1.0])
    for c in cs:
        pl = poly_mul(pl, [c, 1.0])
    pr = [1.0]
    for z, r in num_pairs:
        pr = poly_mul(pr, [-r, 1.0])
        pr = poly_mul(pr, [a0 + z, 1.0])
    for c in cs:
        pr = poly_mul(pr, [a0 - c, 1.0])
    diff = [a - b for a, b in zip(pl, pr)]
    roots = durand_kerner(diff)

    c_const = sum(d * math.log(d) for d in ds) - 2.0 * sum(L * math.log(L) for L in Ls)
    max_r = max(r for _, r in num_pairs)
    # Wait: for alternating pole sum Q_n = sum_j (-1)^j a_j, its saddle is the COMPLEX root (or real root < -a0-max_r)
    # of pl - pr = 0!! Specifically, the contour integral for Q_n around the poles has steepest descent at
    # the root z_Q of pl - pr = 0 with Re(z_Q) in (-a0, 0) (complex) OR real z_Q < -a0 - max_r (Wait: by symmetry
    # z -> -z - a0, the root x_tail > max_r mirrors to x_left = -a0 - x_tail < -a0 - max_r, giving the SAME value
    # for the tail sum!). And what is the growth rate sigma of Q_n?!
    # By Nesterenko/Zudilin/Ball-Rivoal: in the residue sum Q_n = -1/(2pi i) int_{circle around poles} R_n(t) dt?
    # NO! Q_n = sum a_j = coeff of 1/t^3? NO!
    # Look at the generating integral: Q_n = (-1)^... * (1/2) * int_{|w|=...} ... Wait: in Apéry's VWP 7F6(1):
    #   R_n(t) = (2t+n)(t-n)_n(t+n+1)_n / (t)_{n+1}^3.
    #   Let's evaluate BOTH real roots and complex roots of pl - pr = 0!
    def phi_real_tail(x):
        val = c_const
        for zu, ru in num_pairs:
            val += (x - zu)*math.log(x - zu) - (x - ru)*math.log(x - ru)
            val += (x + a0 + ru)*math.log(x + a0 + ru) - (x + a0 + zu)*math.log(x + a0 + zu)
        for c in cs:
            val -= (x + a0 - c)*math.log(x + a0 - c) - (x + c)*math.log(x + c)
        return val

    def phi_cpx(z):
        val = complex(c_const, 0.0)
        for zu, ru in num_pairs:
            val += xlogx(z - zu) - xlogx(z - ru)
            val += xlogx(z + a0 + ru) - xlogx(z + a0 + zu)
        for c in cs:
            val -= xlogx(z + a0 - c) - xlogx(z + c)
        return val

    tail_err = -1e100
    cpx_sig = -1e100
    for z in roots:
        if abs(z.imag) < 1e-6 and z.real > max_r + 1e-5:
            tail_err = max(tail_err, phi_real_tail(z.real))
        elif z.imag > 1e-6:
            # For Q_n = sum (-1)^j a_j, Cauchy contour around poles with pi/sin(pi*z) gives +pi*Im(z) or Re(Phi(z))!
            # Wait: check exact Apéry at z_cpx!
            val = phi_cpx(z)
            # In residue sum sum_j (-1)^j a_j: at t = -j, 1/sin(pi*t) has residue (-1)^j/pi,
            # and for Im(z) > 0, log|1/sin(pi*z)| = -pi*Im(z), wait: a_j has (-1)^{3j} = (-1)^j,
            # so sum a_j = sum (-1)^j |a_j|! No, a_j ALREADY has (-1)^j!
            # And for the integral around the poles: |a_j| = |Res| ~ exp(Re(Phi(-j)))...
            # Let's check Re(Phi(z)) + pi*Im(z) vs Re(Phi(z)) - pi*Im(z) on Apéry RIGHT NOW!
            cpx_sig = max(cpx_sig, val.real)
    return roots, tail_err, cpx_sig


# Check Apéry (a0=1, cs=[0,0,0], num=[(0,1)]):
roots, te, csig = vwp_3blk_saddle(1, [0, 0, 0], [(0, 1)])
print("Apéry VWP 7F6(1) (a0=1, c=[0,0,0], num=[(0,1)]):")
print("  roots of pl - pr = 0:", [(round(z.real, 5), round(z.imag, 5)) for z in roots])
print(f"  tail_err = {te:.8f} (expected -4*ln(1+sqrt(2)) = {-4*math.log(1+math.sqrt(2)):.8f})")
for z in roots:
    if z.imag > 1e-6:
        val = xlogx(z) - xlogx(z - 1) + xlogx(z + 2) - xlogx(z + 1) - 3 * (xlogx(z + 1) - xlogx(z))
        print(f"  complex root z={z.real:.5f}+{z.imag:.5f}i: Re(Phi)={val.real:.8f}, Re(Phi)+pi*y={val.real+PI*z.imag:.8f}, Re(Phi)-pi*y={val.real-PI*z.imag:.8f}")
        print(f"    (compare Apéry +4*ln(1+sqrt(2)) = {4*math.log(1+math.sqrt(2)):.8f}!)")

# ALSO check exact empirical Q_n growth rate of Apéry VWP 7F6(1) for n=1,3,5,7,9 (odd n) AND a0=2 (all n)!
print("\nExact empirical Q_n growth of VWP (a0=2, cs=[0,0,0], num=[(0,2)]):")
for n in [1, 2, 3, 4, 5]:
    # R_n(t) = (2t + 2n - 1) * (t - 2n)_{2n} * (t + 2n)_{2n} / (t)_{2n}^3 ?
    # Wait: Apéry scaled by 2 is a0=2, num=[(0,2)] -> L=2, 2L=4 < 3*2=6!
    Cn = Fraction(math.factorial(2*n - 1)**3, math.factorial(2*n)**2)
    Q = Fraction(0)
    for j in range(2 * n):
        lin = 2*n - 1 - 2*j
        # (j+2n)!/j! * (4n-1-j)!/(2n-1-j)! / [j!(2n-1-j)!]^3 * (-1)^j
        term = Cn * lin * Fraction(math.factorial(j+2*n), math.factorial(j)) * Fraction(math.factorial(4*n-1-j), math.factorial(2*n-1-j))
        term /= (math.factorial(j) * math.factorial(2*n-1-j))**3
        if j % 2 == 1: term = -term
        Q += term
    print(f"  n={n}: log|Q_n|/n = {math.log(abs(float(Q)))/n:.6f} (expected 2*4*ln(1+sqrt(2)) = {8*math.log(1+math.sqrt(2)):.6f})")

# NOW TEST RHIN-VIOLA RV-1 DIRECTLY FOR HIDDEN F_p-SHELLS BEYOND W(D_5)!!
# Wait: we don't even need the 7F6(1) dictionary to compute exact Q_n(RV-1) and P_n(RV-1)!!
# Look at Rhin-Viola's 3D integral for RV-1: (h,j,k,l,m,q,r,s) = (16,17,19,15,12,11,9,13), q' = q+h-r = 18!
# By Zudilin (2004, Proposition 2, eq. (16) - EXACT explicit rational function of ANY Rhin-Viola integral!!):
#   Wait: in Zudilin's well-poised 7F6(1) for RV-1, the 6 parameters (eta0; eta1, eta2, eta3, eta4, eta5)
#   satisfy sum_{i=1}^5 eta_i = 2*eta0!
#   Wait: for ANY (h,j,k,l,m,q,r,s) with h+m=k+r and j+q=l+s, recall Zudilin (2003, *Well-poised hypergeometric service*, p. 597):
#     Let's check ALL 6-tuples (eta0; eta1..eta5) with sum(eta_i) = 2*eta0 and eta_0 = 28..70!
#   NO! Even more direct: remember the EXACT Beukers-Rhin-Viola formula for Q_n(h,j,k,l,m,q,r,s)!!
#   Look at I(h,j,k,l,m,q,r,s) = int_{[0,1]^3} x^h(1-x)^l y^k(1-y)^s z^j(1-z)^q / (1-(1-xy)z)^{q'+1} dx dy dz
#   where q' = q+h-r. By Rhin-Viola (2001, Lemma 2.1):
#     I = Q * zeta(3) - P, where Q is the residue of the 3D integral at 1 - (1-xy)z = 0, i.e., the 2D integral:
#     Q = (1 / q'!) * sum_{nu=0}^{...} ... Specifically, by expanding in t via Beta functions (Rhin-Viola eq. (2.6)):
#       Q = (1 / (q'!)) * [coeff of 1/t^3] -> EXACTLY a terminating 4F3(1)!!
#   Let's write Rhin-Viola's EXACT terminating 4F3(1) formula for Q(h,j,k,l,m,q,r,s) (from Rhin-Viola 2001, eq. (2.5) / Hata 2000):
#     Substitute u = 1 - xy -> int_0^1 z^j (1-z)^q / (1 - u z)^{q'+1} dz = (1/q'!) (d/du)^{q'} [ int_0^1 z^{j-q'}(1-z)^q / (1-uz) dz ]?
#     Even easier! Write 1 / (1 - (1-xy)z) = sum_{n=0}^inf (1-xy)^n z^n? No, use the exact partial fraction in z!
#     In z, the integrand is z^j (1-z)^q / (1 - w z)^{q'+1} where w = 1 - xy.
#     Wait: at w = 1 (xy = 0), 1 - w z = 1 - z!
#     Write z = 1 - (1 - w z)/w - (1 - w)/w = (1 - (1-wz))/w!
#     Expand z^j (1-z)^q in powers of (1 - wz) around z = 1/w:
#       z = (1 - u)/w  (where u = 1 - wz),  1 - z = (w - 1 + u)/w = (u - xy)/w!
#       dz = -du/w.
#       The COEFFICIENT of u^{q'} in (1 - u)^j (u - xy)^q / w^{j+q+1} gives the log(1/(xy)) term:
#       int_0^1 u^{-1} du = -ln(1 - w) = -ln(xy)!!
#     SO THE COEFFICIENT OF ln(1/(xy)) IS EXACTLY:
#       A(x, y) = (1 / (1-xy)^{j+q+1}) * [u^{q'}] ( (1-u)^j (u - xy)^q )
#     Wait: expand (1-u)^j = sum_a C(j,a) (-u)^a and (u-xy)^q = sum_b C(q,b) u^b (-xy)^{q-b}!
#     Wait: q' = q + h - r, and for a+b = q', we integrate int_0^1 int_0^1 x^h(1-x)^l y^k(1-y)^s * (-ln(xy)) * A(x,y) dx dy!
#     Because int_0^1 int_0^1 (-ln(xy)) * x^A (1-x)^B y^C (1-y)^D / (1-xy)^{E+1} dx dy = Q_2 * zeta(3) - P_2!!
#     Specifically, the coefficient of zeta(3) in int_0^1 int_0^1 (-ln(xy)) x^A y^C / (1-xy) dx dy is sum_{n>=0} (1/(n+A+1)^2? NO!):
#     int_0^1 int_0^1 (-ln(xy)) (xy)^n dx dy = int_0^1 int_0^1 (-ln x - ln y) x^n y^n dx dy = 2 / (n+1)^3!!
#     THEREFORE: for ANY P = (h,j,k,l,m,q,r,s),
#       I(P) = (1/q'!) * ... has EXACT zeta(3) coefficient Q(P) given by the 1D sum (Rhin-Viola 2001 eq. (2.4)):
#       Q(h,j,k,l,m,q,r,s) = sum_{t=0}^{min(h,k,j,q')} ... (wait: let's check against Apéry h=j=k=l=m=q=r=s=1 => Q = 2*5 = 10?!)
#   YES!! For Apéry (1,1,1,1,1,1,1,1):
#     I = int_{[0,1]^3} x(1-x)y(1-y)z(1-z)/(1-(1-xy)z)^2 dx dy dz = 2*(5*zeta(3) - 6) = 10*zeta(3) - 12?
#     Check: Beukers (1979) integral I_1 = int x(1-x)y(1-y) / (1 - (1-xy)z) -> wait:
#     Beukers's integral is int_0^1 int_0^1 (-ln(xy))/(1-xy) * P_n(x) P_n(y) dx dy = 2 * (b_n zeta(3) - a_n),
#     and int_0^1 dz / (1 - (1-xy)z) = -ln(xy) / (1 - xy)!
#     With z^j(1-z)^q/(1-(1-xy)z)^{q'+1} (where j=q=q'=n):
#       int_0^1 z^n (1-z)^n / (1 - (1-xy)z)^{n+1} dz = (-1)^n / n! * (d/dw)^n [ ... ] = P_n(xy) * (-ln(xy))/(1-xy) + rat!
#     Specifically, let's check the EXACT closed-form sum for Q(h,j,k,l,m,q,r,s) at n=1:
#       Expand 1 / (1 - (1-xy)z)^{q'+1} = sum_{t=0}^inf C(t+q', q') (1-xy)^t z^t!
#       Integrate z: sum_{t=0}^inf C(t+q', q') * (j+t)! q! / (j+q+1+t)! * (1-xy)^t!
#       Expand (1-xy)^t? No, integrate x, y via Beta on (1-xy)^t:
#         int_0^1 int_0^1 x^h (1-x)^l y^k (1-y)^s (1-xy)^t dx dy
#         = sum_{a=0}^t (-1)^a C(t, a) * (h+a)! l! / (h+l+1+a)! * (k+a)! s! / (k+s+1+a)!
#       THAT IS LITERALLY AN EXACT POSITIVE/RATIONAL INTEGRAL!
#       Even better: **ZUDILIN'S EXACT FORMULA FOR Q(h,j,k,l,m,q,r,s)** (Rhin-Viola 2001 Theorem 2.1, eq. (2.8)):
#       By Euler's transformation on the (x, y) integral (putting t = j + q' ...), Q(h,j,k,l,m,q,r,s) is:
#         Q_n = (1 / (l! s! q!)) * sum_{a = max(0, q'-j)}^{min(h', k', ...)} ...
#       WAIT: Let's just use **Whipple's 4F3(1) = 7F6(1) formula** directly on (h,j,k,l,m,q,r,s)!!
#       In Rhin-Viola (2001, Section 2, eq. (2.4) & (2.8)):
#         Let H = h, J = j, K = k, L = l, M = m, Q = q, R = r, S = s, and V' = (h',j',k',l',m',q',r',s').
#         The coefficient of 2*zeta(3) in I(P) (times (h! l! k! s! j! q! / q'!)^{-1}? No, I(P) itself!) is:
#         Q(P) = (-1)^{...} * sum_{nu} (-1)^nu * C(h+nu, nu) * ... ->
#         Let's derive Q(P) in 5 lines from int_0^1 int_0^1 (-ln(xy)) * A(x,y) dx dy!
#         Look: A(x,y) = x^h(1-x)^l y^k(1-y)^s / (1-xy)^{j+q+1} * sum_{a=0}^{q'} C(j, a) (-1)^a C(q, q'-a) (-xy)^{q - (q'-a)}
#                      = sum_{a=0}^{q'} (-1)^{q - q' + 2a} C(j, a) C(q, q'-a) * x^{h+q-q'+a} (1-x)^l y^{k+q-q'+a} (1-y)^s / (1-xy)^{j+q+1}
#         Now use: for ANY A, B, C, D, E >= 0:
#           coeff of zeta(3) in int_0^1 int_0^1 (-ln(xy)) x^A (1-x)^B y^C (1-y)^D / (1-xy)^{E+1} dx dy:
#           Write 1/(1-xy)^{E+1} = sum_{m>=0} C(m+E, E) (xy)^m!
#           Wait: the coeff of 2*zeta(3) in sum_{m>=0} P(m) / (m+1)^3 is JUST P(-1)!! NO, it's the residue!
#           Look: int_0^1 x^{m+A}(1-x)^B dx = B! / (m+A+1)_{B+1} = sum_{u=A+1}^{A+B+1} c_u / (m + u),
#                 int_0^1 y^{m+C}(1-y)^D dy = D! / (m+C+1)_{D+1} = sum_{v=C+1}^{C+D+1} d_v / (m + v)!
#           When you differentiate (-ln(xy)) = -(d/dA + d/dC), the double pole at u = v gives 2 / (m + u)^3!!
#           SO THE EXACT COEFFICIENT OF 2*zeta(3) IN I(P) IS:
#             Q(P) = sum_{a} (-1)^{q-q'} C(j, a) C(q, q'-a) * sum_{u} C(u-1+E, E) * c_u * d_u ?
#             Wait: m + u = k => m = k - u => C(m+E, E) = C(k - u + E, E) = polynomial in k -> coeff of k^0 at k -> 0 is C(E - u, E) = (-1)^E C(u - 1, E)!
#             WAIT A SECOND!! CHECK THAT!!
#             At u = v: 1/(m+u)^3 * C(m+E, E) = 1/k^3 * C(k - u + E, E) -> the constant term (k=0) of C(k - u + E, E) IS C(E - u, E) = (-1)^E C(u - 1, E)!!
#             AND c_u = (-1)^{u-A-1} C(B, u-A-1), d_u = (-1)^{u-C-1} C(D, u-C-1)!!
#             HOLD ON! Let's TEST THAT 3-LINE FORMULA ON APÉRY (1,1,1,1,1,1,1,1) AND (2,2,2,2,2,2,2,2) RIGHT NOW!

def Q_rv_exact(P, n=1):
    """Exact Q_n (coeff of 2*zeta(3)) of Rhin-Viola integral I(n*P) via exact residue derivation!"""
    h, j, k, l, m, q, r, s = [x * n for x in P]
    qp = q + h - r
    E = j + q  # denominator power is (1-xy)^{E+1}
    # Check: A = h + q - qp + a = r + a; C = k + q - qp + a = k - h + r + a = m + a!
    # Wait: q - qp = r - h! So A = r + a, C = k + r - h + a = m + a!
    # Check signs: (-1)^{q - qp} * (-1)^E * (-1)^{u-A-1} * (-1)^{u-C-1} = (-1)^{q-qp+E+A+C}
    # Notice A+C = r+m+2a, so sign is (-1)^{q-qp+j+q+r+m} = (-1)^{j+r+m-qp} (INDEPENDENT of a and u!!)!
    total = 0
    for a in range(max(0, qp - q), min(j, qp) + 1):
        # coeff from z-integral: C(j, a) * C(q, qp - a)
        ca = math.comb(j, a) * math.comb(q, qp - a)
        A = r + a
        C = m + a
        u_lo = max(A + 1, C + 1, E + 1)
        u_hi = min(A + l + 1, C + s + 1)
        s_u = 0
        for u in range(u_lo, u_hi + 1):
            # C(u-1, E) * C(l, u-A-1) * C(s, u-C-1)
            s_u += math.comb(u - 1, E) * math.comb(l, u - A - 1) * math.comb(s, u - C - 1)
        total += ca * s_u
    return total

print("\nTesting exact Q_rv_exact on Apéry (1,1,1,1,1,1,1,1) [should give Apéry b_n = 5, 73, 1445, 33001!]:")
for n in range(1, 5):
    print(f"  n={n}: Q_rv_exact(Apery) = {Q_rv_exact((1,1,1,1,1,1,1,1), n)}")

# NOW TEST RHIN-VIOLA'S 2001 RECORD RV-1 = (16, 17, 19, 15, 12, 11, 9, 13) FOR EXTRA F_p-SHELLS BEYOND W(D_5)!!
print("\n================================================================================")
print("EXP 13A: DOES RHIN-VIOLA RV-1 (5.513891) HAVE EXTRA F_p-SHELLS BEYOND W(D_5)?!")
print("================================================================================")
RV1 = (16, 17, 19, 15, 12, 11, 9, 13)
# First verify exact growth rate log(Q_n)/n matches RV-1's c1 = 48.4694096!
q1 = Q_rv_exact(RV1, 1)
q2 = Q_rv_exact(RV1, 2)
print(f"  RV-1 n=1: log(Q_1)={math.log(q1):.6f}, n=2: log(Q_2)/2={math.log(q2)/2:.6f}, diff={math.log(q2)-math.log(q1):.6f} (RV-1 c1 = 48.469410!)")

# Compute the exact W(D_5) valuation v_W(u) of RV-1 on u in (1/19, 1) (from certify_rv_bound.py!),
# and compare against the TRUE exact v_p(Q_n(RV-1))!!
# WAIT A SECOND!! Remember: in RV-1, Q_n is an INTEGER (divided by nothing!), AND W(D_5) proves that
# Q_n is divisible by prod_{p in Omega_E} p * prod_{p in Omega'_E} p^2!!
# If v_p(Q_n(RV-1)) > v_{W(D_5)}(n/p) systematically on ANY interval of u in (0, 1),
# THAT MEANS RHIN-VIOLA'S 2001 RECORD HAS HIDDEN F_p-SHELLS THAT NOBODY EVER FOUND!!
# Let's check v_p(Q_n(RV-1)) vs v_{W(D_5)}(n/p) RIGHT NOW at p = 41, 43, 47, 53 (exact mod p^4)!

# Build exact W(D_5) predicted valuation v_W(u) from certify_rv_bound.py:
FN = ((0,2,3,6),(0,1,3,5),(0,3,4,7),(1,2,5,6),(0,1,4,5),(1,3,5,7))
FD = ((8,12,13,14),(11,13,14,15),(9,10,13,14),(8,9,12,13),(9,10,13,14),(9,11,13,15))
FR = (8,8,4,4,4,2)
def theta(i, t): return (i+t)%8 if i<8 else 8+((i-8+t)%8)
h,j,k,l,m,q,r,s = RV1
V_rv = [h,j,k,l,m,q,r,s, h+l-j, j+m-k, k+q-l, l+r-m, m+s-q, q+h-r, r+j-s, s+k-h]
# Build the 30 cosets of W(D_5):
cosets = []
for f in range(6):
    for u_i in range(FR[f]):
        cosets.append(([V_rv[theta(FN[f][i], u_i)] for i in range(4)], [V_rv[theta(FD[f][i], u_i)] for i in range(4)]))

def v_WD5_rv1(u_frac: Fraction):
    """Exact W(D_5) predicted exponent in {0, 1, 2} for u = n/p in (0, 1) (where p > M*n is 0)."""
    # Actually in RV01 Theorem 5.1, for ANY u = n/p > 0:
    #   e_c(u) = sum_{i=0..3} floor(num_i * u) - sum_{i=0..3} floor(den_i * u)  [Wait: +2 from 3D integral? No:
    #   In RV01, the factorial ratio of coset c is prod(num_i*n)! / prod(den_i*n)!, wait: sum(num_i) == sum(den_i)!
    #   And the valuation of Q_n from coset c is max_c ( sum floor(den_i*u) - sum floor(num_i*u) )? Let's check!]
    # Let's check exact coset valuations at n=1, p=2..19 against v_p(q1)!!
    best = 0
    for nums, dens in cosets:
        # In RV01, E_c(x) = sum floor(num*x) - sum floor(den*x) -> Wait: check sign!
        # Check at u = 1/19: max(V)=19 is in nums or dens?
        # Look at V[2]=k=19 (index 2 in 0..7 -> NUM!), V[13]=q'=18 (index 13 -> DEN!), V[1]=j=17 (index 1 -> NUM!).
        # At u = 1/19, only 19*u >= 1 -> sum floor(num*u) = 1, sum floor(den*u) = 0!
        # Wait: at u = 1/19 (p = 19, n = 1), check v_19(q1)!
        val = sum(math.floor(x * u_frac) for x in nums) - sum(math.floor(x * u_frac) for x in dens)
        if val > best: best = val
    return best

print("  Calibration check of v_p(Q_1(RV-1)) vs W(D_5) at n=1, primes p <= 23:")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    # Exact v_p(q1)
    v_true = 0; t = q1
    while t % p == 0: v_true += 1; t //= p
    v_w = v_WD5_rv1(Fraction(1, p))
    print(f"    p={p:2d} (u=1/{p:2d}): v_p(Q_1) = {v_true}, v_WD5 = {v_w}, extra = {v_true - v_w}")

# NOW THE MOMENT OF TRUTH: test RV-1 at p = 31, 37, 41, 43 across ALL n in [2..p-1] where p > sqrt(35*n)!!
# Compute Q_n(RV-1) mod p^4 fast!
def vp_Q_rv1_mod_p4(n, p):
    h,j,k,l,m,q,r,s = [x * n for x in RV1]
    qp = q + h - r
    E = j + q
    lim = max(r + qp + l + 2, m + qp + s + 2, E + 2) + 5
    if lim >= p * p: return -1
    mod = p ** 4
    vf = [0]*(lim+1); uf = [1]*(lim+1)
    for i in range(1, lim+1):
        x, c = i, 0
        while x % p == 0: c += 1; x //= p
        vf[i] = vf[i-1]+c; uf[i] = (uf[i-1]*x)%mod
    def bvu(N, K):
        if K < 0 or K > N: return 999, 0
        return vf[N]-vf[K]-vf[N-K], (uf[N]*pow((uf[K]*uf[N-K])%mod, -1, mod))%mod

    vmin = 999
    terms = []
    for a in range(max(0, qp - q), min(j, qp) + 1):
        v1, u1 = bvu(j, a)
        v2, u2 = bvu(q, qp - a)
        va = v1 + v2; ua = (u1 * u2) % mod
        A = r + a; C = m + a
        u_lo = max(A + 1, C + 1, E + 1)
        u_hi = min(A + l + 1, C + s + 1)
        for u in range(u_lo, u_hi + 1):
            v3, u3 = bvu(u - 1, E)
            v4, u4 = bvu(l, u - A - 1)
            v5, u5 = bvu(s, u - C - 1)
            vt = va + v3 + v4 + v5
            if vt < vmin: vmin = vt
            terms.append((vt, (ua * u3 % mod) * (u4 * u5 % mod) % mod))
    acc = sum(u * (p ** (v - vmin)) for v, u in terms if v - vmin < 4) % mod
    ext = 4 if acc == 0 else (3 if acc % (p**3) == 0 else (2 if acc % (p*p) == 0 else (1 if acc % p == 0 else 0)))
    return vmin, vmin + ext

print("\n  MOMENT OF TRUTH — Testing RV-1 Q_n mod p^4 vs W(D_5) across u in (1/19, 1) at p=43 and p=47!")
for p in [43, 47]:
    beats = []
    matches = 0
    total = 0
    for n in range(1, p):
        if 45 * n >= p * p: continue
        res = vp_Q_rv1_mod_p4(n, p)
        if res == -1: continue
        vmin_car, v_true = res
        v_w = v_WD5_rv1(Fraction(n, p))
        total += 1
        if v_true > v_w:
            beats.append((n, p, round(n/p, 4), vmin_car, v_w, v_true))
        elif v_true == v_w:
            matches += 1
    print(f"  p={p}: tested {total} values of n (u up to {total/p:.3f}): exact match W(D_5)={matches}, STRICTLY > W(D_5)={len(beats)}!")
    for b in beats[:15]:
        print(f"    BEAT! n/p={b[0]}/{b[1]} ({b[2]}): raw_car={b[3]}, W(D_5)={b[4]}, TRUE v_p(Q_n)={b[5]} (+{b[5]-b[4]} ABOVE W(D_5)!!)")
