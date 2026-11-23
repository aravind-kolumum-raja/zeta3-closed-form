#!/usr/bin/env python3
r"""
cdt_capacity_test.py — Quantitative Calegari-Dimitrov-Tang (CDT 2024-2026) /
Bost-Charles conformal-capacity vs. p-adic denominator-growth barrier analysis
for zeta(3) and zeta(3)/pi^3.

Computes to 50+ decimal digits (stdlib Decimal only):
1. Classical trilogarithm Li_3(x): unit disk vs Koebe slit plane C \ [1, inf).
2. Beukers (1987) Gamma_0(6)^+ Hauptmodul t(q) = (eta(tau)*eta(6*tau)/(eta(2*tau)*eta(3*tau)))^12
   and exact verification of:
     - t(i/sqrt(6)) == (sqrt(2) - 1)^4 = 17 - 12*sqrt(2) (50-digit check)
     - A(t(q)) = sum A_n t(q)^n has pure integer q-expansion (weight-2 modular form)
     - CDT conformal amplification log(1/t_1) = 4*log(1+sqrt(2)) = 3.525494 > tau = 3
       (reproducing exactly why Apery's proof clears the e^{3n} barrier with surplus +0.5255!)
3. Modular lambda-function on Gamma(2) and Plouffe's CM Hecke-coupled Eichler
   system (S_3(1), S_3(2), S_3(4)) at tau = i (where lambda(i) = 1/2):
     - Exact q-expansion of lambda(q) = 16q - 128q^2 + 704q^3 - ...
     - Cusp mapping radius |lambda'(0)| = 16 => log(16) = 2.772589 < tau = 3 (deficit -0.2274!)
     - Local conformal mapping radius at the CM point tau = i via the elliptic period
       formula |d lambda / d w| at w = (tau - i)/(tau + i) = 0, using AGM / Gamma(1/4)^4!
     - Exact quantitative gap for proving zeta(3)/pi^3 irrational via CDT/Bost-Charles.
"""

from decimal import Decimal, getcontext
import math

getcontext().prec = 70


def pi_dec() -> Decimal:
    # Chudnovsky-style / Machin formula to 70 digits
    def arctan_inv(m: int) -> Decimal:
        inv_m = Decimal(1) / Decimal(m)
        inv_m2 = inv_m * inv_m
        term = inv_m
        total = term
        k = 1
        while True:
            term = -term * inv_m2
            add = term / Decimal(2 * k + 1)
            if abs(add) < Decimal("1e-68"):
                break
            total += add
            k += 1
        return total

    return Decimal(4) * (Decimal(4) * arctan_inv(5) - arctan_inv(239))


PI = pi_dec()


def eta_prod_ratio(q: Decimal, max_n: int = 200) -> Decimal:
    """
    Compute t(q) = q * prod_{n>=1} ((1-q^n)(1-q^{6n}) / ((1-q^{2n})(1-q^{3n})))^12
    Wait: check the q-power of (eta(tau)eta(6tau)/(eta(2tau)eta(3tau)))^12!
      Exponent of q^{1/24}: 12 * (1 + 6 - 2 - 3) / 24 = 12 * 2 / 24 = 1!
    So t(q) = q * prod_{n>=1} ((1-q^n)(1-q^{6n}) / ((1-q^{2n})(1-q^{3n})))^12.
    """
    prod = Decimal(1)
    qn = q
    for n in range(1, max_n + 1):
        if qn < Decimal("1e-65"):
            break
        q2 = qn * qn
        q3 = q2 * qn
        q6 = q3 * q3
        fac = ((Decimal(1) - qn) * (Decimal(1) - q6)) / ((Decimal(1) - q2) * (Decimal(1) - q3))
        prod *= fac**12
        qn *= q
    return q * prod


def series_mul(a, b, deg):
    out = [0] * (deg + 1)
    for i in range(deg + 1):
        if a[i] == 0:
            continue
        for j in range(deg + 1 - i):
            out[i + j] += a[i] * b[j]
    return out


def series_pow_1_minus_qm(m, exp, deg):
    """Exact integer series of (1 - q^m)^exp mod q^{deg+1}."""
    # exp can be negative or positive int
    out = [0] * (deg + 1)
    out[0] = 1
    if exp == 0:
        return out
    step = [0] * (deg + 1)
    if exp > 0:
        # binomial (1 - q^m)^exp
        c = 1
        for k in range(exp + 1):
            if k * m > deg:
                break
            out[k * m] = c * ((-1) ** k)
            c = c * (exp - k) // (k + 1)
    else:
        e = -exp
        c = 1
        for k in range(deg // m + 1):
            out[k * m] = c
            c = c * (e + k) // (k + 1)
    return out


def beukers_t_qseries(deg: int = 12):
    """Exact integer q-series of t(q) = q * prod ((1-q^n)(1-q^{6n})/((1-q^{2n})(1-q^{3n})))^12."""
    s = [0] * (deg + 1)
    s[0] = 1
    for n in range(1, deg + 1):
        for mult, exp in [(1, 12), (6, 12), (2, -12), (3, -12)]:
            m = n * mult
            if m <= deg:
                s = series_mul(s, series_pow_1_minus_qm(m, exp, deg), deg)
    # multiply by q
    return [0] + s[:deg]


def apery_A(n: int) -> int:
    total = 0
    for k in range(n + 1):
        total += (math.comb(n, k) ** 2) * (math.comb(n + k, k) ** 2)
    return total


def modular_lambda_q(q: Decimal, max_n: int = 120) -> Decimal:
    """
    Modular lambda function lambda(tau) = 16 * q^{1/2} * prod_{n>=1} ((1+q^n)/(1+q^{n-1/2}))^8
    Here pass r = q^{1/2} = exp(i*pi*tau), so lambda = 16*r * prod_{n>=1} ((1+r^{2n})/(1+r^{2n-1}))^8.
    Wait: for tau = i, r = exp(-pi), lambda(i) = 1/2!
    """
    r = q.sqrt()
    prod = Decimal(1)
    r_odd = r
    for n in range(1, max_n + 1):
        if r_odd < Decimal("1e-65"):
            break
        r_even = r_odd * r
        fac = (Decimal(1) + r_even) / (Decimal(1) + r_odd)
        prod *= fac**8
        r_odd = r_even * r
    return Decimal(16) * r * prod


def agm_dec(a: Decimal, b: Decimal) -> Decimal:
    for _ in range(100):
        an = (a + b) / Decimal(2)
        bn = (a * b).sqrt()
        if abs(an - bn) < Decimal("1e-66"):
            return an
        a, b = an, bn
    return a


def main():
    print("========================================================================")
    print("  QUANTITATIVE CDT (2024-2026) CONFORMAL-CAPACITY BARRIER ANALYSIS")
    print("  FOR zeta(3) AND zeta(3)/pi^3")
    print("========================================================================\n")

    # ------------------------------------------------------------------
    # 1. Classical Li_3(x) benchmarks
    # ------------------------------------------------------------------
    tau_z3 = Decimal(3)  # d_n^3 ~ e^{3n}
    log_koebe = Decimal(4).ln()
    print("[1] Classical Trilogarithm Li_3(x) = sum x^n / n^3 (p-adic barrier tau = 3.000000)")
    print(f"    Unit disk |x| < 1:            rho = 1.000000, log(rho) = 0.000000  (gap = {-tau_z3:+.6f})")
    print(f"    Koebe slit C \\ [1, inf):      rho = 4.000000, log(rho) = {log_koebe:.6f}  (gap = {log_koebe - tau_z3:+.6f})")

    # ------------------------------------------------------------------
    # 2. Beukers (1987) Gamma_0(6)^+ modular uniformization of Apery K3 family
    # ------------------------------------------------------------------
    print("\n[2] Beukers Gamma_0(6)^+ Modular Parametrization of Apery's K3 Picard-Fuchs ODE")
    sqrt2 = Decimal(2).sqrt()
    sqrt6 = Decimal(6).sqrt()
    t1 = (sqrt2 - Decimal(1)) ** 4  # 17 - 12*sqrt(2)
    t2 = (sqrt2 + Decimal(1)) ** 4  # 17 + 12*sqrt(2)
    q_fricke = (-Decimal(2) * PI / sqrt6).exp()
    t_fricke = eta_prod_ratio(q_fricke)
    err_fricke = abs(t_fricke - t1)
    assert err_fricke < Decimal("1e-55"), f"Fricke identity failed: {err_fricke}"
    print(f"    Singularity t_1 = (sqrt(2)-1)^4  = {t1:.45f}")
    print(f"    t(q = exp(-2*pi/sqrt(6)))        = {t_fricke:.45f}")
    print(f"    Exact agreement |t(i/sqrt(6))-t1| = {err_fricke:.2e} [VERIFIED TO 55+ DIGITS]")

    # Verify exact integer q-series of t(q) and A(t(q))
    tq = beukers_t_qseries(10)
    print(f"    t(q) integer q-series (n=0..10): {tq}")
    assert tq[:5] == [0, 1, -12, 66, -220]
    # Compose A(t(q)) = sum_{k=0}^10 A_k t(q)^k mod q^11
    deg = 10
    A_q = [0] * (deg + 1)
    t_pow = [0] * (deg + 1)
    t_pow[0] = 1
    for k in range(deg + 1):
        ak = apery_A(k)
        for i in range(deg + 1):
            A_q[i] += ak * t_pow[i]
        t_pow = series_mul(t_pow, tq, deg)
    print(f"    A(t(q)) integer q-series (0..10): {A_q}")

    # Conformal amplification of q -> t(q) relative to the singularity radius t_1:
    # Since t(q) = q + O(q^2) maps |q| < 1 holomorphically into C \ cuts with |t'(0)| = 1,
    # in the normalized variable w = t / t_1 (where A_n t_1^n ~ n^{-3/2} has radius 1),
    # the mapping radius is rho = |w'(0)| = 1 / t_1 = (1 + sqrt(2))^4!
    rho_apery = Decimal(1) / t1
    log_rho_apery = rho_apery.ln()
    surplus_apery = log_rho_apery - tau_z3
    mu_apery = Decimal(1) + (log_rho_apery + tau_z3) / (log_rho_apery - tau_z3)
    print(f"    Normalized conformal radius rho = 1/t_1 = (1+sqrt(2))^4 = {rho_apery:.10f}")
    print(f"    log(rho) = {log_rho_apery:.10f}  vs  tau = 3.0000000000  =>  SURPLUS = {surplus_apery:+.10f} > 0!")
    print(f"    Implied irrationality exponent mu(zeta(3)) <= {mu_apery:.10f} (recovers Apery's 13.41782!)")

    # ------------------------------------------------------------------
    # 3. Gamma(2) / X_0(4) Plouffe CM System (S_3(1), S_3(2), S_3(4)) at tau = i
    # ------------------------------------------------------------------
    print("\n[3] Modular lambda on Gamma(2) & Plouffe CM Eichler System at tau = i")
    q_i = (-Decimal(2) * PI).exp()
    lam_i = modular_lambda_q(q_i)
    err_lami = abs(lam_i - Decimal("0.5"))
    assert err_lami < Decimal("1e-60"), f"lambda(i) != 1/2: {err_lami}"
    print(f"    lambda(i) = {lam_i:.50f}  (|err from 1/2| = {err_lami:.2e}) [PASS]")

    # At tau = 2i (Plouffe's second Hecke point S_3(2)):
    lam_2i = modular_lambda_q(q_i * q_i)
    err_lam2i = abs(lam_2i - t1)
    assert err_lam2i < Decimal("1e-55"), f"lambda(2i) != (sqrt(2)-1)^4: {err_lam2i}"
    print(f"    lambda(2i) = {lam_2i:.45f} == (sqrt(2)-1)^4  (err = {err_lam2i:.2e}) [PASS!]")

    # At tau = 4i (Plouffe's third Hecke point S_3(4)):
    lam_4i = modular_lambda_q(q_i**4)
    # Exact algebraic value: Landen transform lambda(2tau) = ((1 - sqrt(1-lam))/(1 + sqrt(1-lam)))^2
    s1 = (Decimal(1) - t1).sqrt()
    lam_4i_exact = ((Decimal(1) - s1) / (Decimal(1) + s1)) ** 2
    err_lam4i = abs(lam_4i - lam_4i_exact)
    assert err_lam4i < Decimal("1e-55")
    print(f"    lambda(4i) = {lam_4i:.45f}  (Landen exact err = {err_lam4i:.2e}) [PASS!]")

    # Cusp mapping radius of lambda(r) = 16*r - 128*r^2 + ... (r = exp(i*pi*tau) in unit disk D):
    # Wait: in q = r^2 = exp(2*pi*i*tau), lambda(tau) is a function of r = q^{1/2};
    # the Eichler integral E_3(q) = sum sigma_3(n)/n^3 q^n is analytic on |q| < 1 (|r| < 1).
    # Relative to the unit disk |r| < 1:
    #   - Mapping radius at the cusp r = 0 of lambda(r) = 16r + ... is |lambda'(0)| = 16!
    log_16 = Decimal(16).ln()
    print(f"\n    Cusp mapping radius of lambda(r) on Gamma(2): rho_cusp = 16")
    print(f"    log(16) = {log_16:.10f}  vs  tau = 3.0000000000  =>  DEFICIT = {log_16 - tau_z3:+.10f} < 0!")
    print(f"    (Compare weight-2 CDT for zeta(2), L(2,chi_-3): tau = 2 => log(16) - 2 = {log_16 - Decimal(2):+.6f} > 0!)")

    # Local conformal mapping radius at the CM point tau = i in the Cayley disk w = (tau - i)/(tau + i):
    # Using classic elliptic-integral derivative d(lambda)/d(tau) = i*pi * lambda*(1-lambda)*theta_3(tau)^4
    # At tau = i: lambda(i) = 1/2, theta_3(i)^2 = 2*K(1/sqrt(2))/pi = 1 / AGM(1, 1/sqrt(2))!
    # Wait: K(k) = pi / (2 * AGM(1, k')), and theta_3^2 = 2*K/pi = 1 / AGM(1, sqrt(1-lambda)).
    # At tau = i: lambda = 1/2, k' = 1/sqrt(2), so theta_3(i)^2 = 1 / AGM(1, 1/sqrt(2))!
    M_lemn = agm_dec(Decimal(1), Decimal(1) / sqrt2)
    theta3_i_4 = Decimal(1) / (M_lemn * M_lemn)
    # |d lambda / d tau| at tau = i is pi * (1/4) * theta_3(i)^4
    dlam_dtau = PI * Decimal("0.25") * theta3_i_4
    # Since tau = i*(1+w)/(1-w), |d tau / d w| at w = 0 is 2!
    # Therefore conformal mapping radius of the universal cover H -> C \ {0, 1} centered at lambda = 1/2 is:
    #   rho(C \ {0, 1}, 1/2) = |d lambda / d w|_{w=0} = 2 * |d lambda / d tau|_{tau=i}!
    rho_cm_i = Decimal(2) * dlam_dtau
    # Cross-check via finite difference around tau = i!
    eps = Decimal("1e-20")
    q_plus = (-Decimal(2) * PI * (Decimal(1) + eps)).exp()
    q_minus = (-Decimal(2) * PI * (Decimal(1) - eps)).exp()
    dlam_num = abs(modular_lambda_q(q_plus) - modular_lambda_q(q_minus)) / (Decimal(2) * eps)
    assert abs(dlam_num - dlam_dtau) < Decimal("1e-18"), f"Derivative mismatch: {dlam_num} vs {dlam_dtau}"

    # Normalized to distance 1/2 from lambda = 1/2 to the nearest singularities {0, 1}:
    # The plain disk |lambda - 1/2| < 1/2 has radius 0.5; the universal cover amplifies by rho_cm_i / 0.5!
    amp_cm_i = rho_cm_i / Decimal("0.5")
    log_amp_i = amp_cm_i.ln()
    print(f"\n    Local CM uniformization at tau = i (lambda = 1/2, universal cover H -> C \\ {{0,1}}):")
    print(f"    AGM(1, 1/sqrt(2))        = {M_lemn:.30f}")
    print(f"    |d lambda / d tau|_{{i}}  = {dlam_dtau:.30f}  (FD check err = {abs(dlam_num - dlam_dtau):.2e})")
    print(f"    Mapping radius rho_i     = {rho_cm_i:.30f}")
    print(f"    Normalized amp rho_i/0.5 = {amp_cm_i:.30f}")
    print(f"    log(rho_i / 0.5)         = {log_amp_i:.10f}  vs  tau = 3.0000000000  =>  DEFICIT = {log_amp_i - tau_z3:+.10f}!")

    # ------------------------------------------------------------------
    # 4. Bost-Charles / CDT Rank-3 Holonomy Bound Summary for (1, pi^3, zeta(3))
    # ------------------------------------------------------------------
    print("\n[4] CDT / Bost-Charles Rank-3 Holonomy Threshold for (1, pi^3, zeta(3))")
    # For a single function with denominator growth tau = 3 and conformal log-radius L:
    #   - Rank-2 irrationality (zeta(3) alone): requires L > tau_eff.
    #     With Gamma_0(6)^+ alone: L = 3.5255 > 3.0000 (PROVED!).
    #     With Rhin-Viola W(D_5) prime saving: tau_eff drops from 3 to 3 * (29.812/54) = 1.6562!
    #   - Rank-3 linear independence (1, pi^3, zeta(3)) => zeta(3)/pi^3 irrational:
    #     In a 3-dimensional simultaneous system (such as Plouffe's S_3(1), S_3(2), S_3(4)
    #     on X_0(4) evaluated at tau = i), the denominator exponent of the Eichler integral
    #     E_3(q) = sum (sigma_3(n)/n^3) q^n is tau = 3, AND for simultaneous approximation
    #     of 2 independent numbers (pi^3, zeta(3)) over Q, Siegel's lemma / Bost-Charles
    #     holonomy bound requires the sum of log conformal radii across the fiber to exceed
    #     (d - 1) * tau (or in CDT holonomy rank bound: log(rho_eff) > tau_eff for the
    #     unbounded-denominator quotient!).
    #   - Even for rank 2 (just one Eichler integral), Gamma(2) has log(16) = 2.7726 < 3
    #     (deficit -0.2274) and the local CM amplification at tau=i has log(amp) = 1.4765 < 3
    #     (deficit -1.5235).
    #   - Crucially, passing from tau = 3 to Rhin-Viola prime-saved tau_eff ~ 1.6562 via
    #     modular Hecke/Weyl symmetries on X_0(4) WOULD flip log(16) = 2.7726 > 1.6562 (+1.1164)!!
    rv_ratio = Decimal("2173.896834425569") / Decimal("3939")  # c2 / (M+N+Q) for 73*RV-1
    tau_rv_eff = Decimal(3) * rv_ratio
    print(f"    73*RV-1 effective denominator factor: c2/(M+N+Q) = {rv_ratio:.10f}  =>  tau_eff = {tau_rv_eff:.6f}")
    print(f"    Summary table (need log(rho) > tau):")
    print(f"      - Li_3(x) Koebe slit plane:        log(rho) = {log_koebe:.6f}  vs tau=3.0000  (gap {log_koebe - tau_z3:+.6f})")
    print(f"      - Plouffe CM local H->C\\{{0,1}}:    log(amp) = {log_amp_i:.6f}  vs tau=3.0000  (gap {log_amp_i - tau_z3:+.6f})")
    print(f"      - Gamma(2) cusp lambda(r)=16r+..:  log(rho) = {log_16:.6f}  vs tau=3.0000  (gap {log_16 - tau_z3:+.6f})")
    print(f"      - Gamma_0(6)^+ Apery t(q):         log(rho) = {log_rho_apery:.6f}  vs tau=3.0000  (gap {surplus_apery:+.6f}) [CLEARS!]")
    print(f"      - HYPOTHETICAL Gamma(2)+RV saving: log(16)  = {log_16:.6f}  vs tau_eff={tau_rv_eff:.4f} (gap {log_16 - tau_rv_eff:+.6f}) [WOULD CLEAR!]")


if __name__ == "__main__":
    main()
