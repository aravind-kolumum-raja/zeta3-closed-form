use std::f64::consts::PI;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Spec {
    rs: Vec<i32>,
    ls: Vec<i32>,
    s1: i32,
    d1: i32,
    s2: i32,
    d2: i32,
}

#[derive(Clone, Copy, Debug)]
struct Cpx {
    re: f64,
    im: f64,
}

impl Cpx {
    fn new(re: f64, im: f64) -> Self { Self { re, im } }
    fn add(self, o: Self) -> Self { Self::new(self.re + o.re, self.im + o.im) }
    fn sub(self, o: Self) -> Self { Self::new(self.re - o.re, self.im - o.im) }
    fn mul(self, o: Self) -> Self { Self::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re) }
    fn div(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Self::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
    fn abs(self) -> f64 { self.re.hypot(self.im) }
    fn ln(self) -> Self { Self::new(self.abs().ln(), self.im.atan2(self.re)) }
    fn xlogx(self) -> Self {
        if self.abs() < 1e-14 { Self::new(0.0, 0.0) } else { self.mul(self.ln()) }
    }
}

fn poly_mul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut r = vec![0.0; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] += x * y;
        }
    }
    r
}

fn durand_kerner(mut coeffs: Vec<f64>) -> Vec<Cpx> {
    while coeffs.len() > 1 && coeffs.last().unwrap().abs() < 1e-10 {
        coeffs.pop();
    }
    let d = coeffs.len() - 1;
    if d == 0 { return vec![]; }
    let lead = coeffs[d];
    let c: Vec<f64> = coeffs.iter().map(|&x| x / lead).collect();
    let max_c = c[..d].iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let rad = (1.0 + max_c).min(60.0);
    let mut roots: Vec<Cpx> = (0..d)
        .map(|k| {
            let ang = 2.0 * PI * (k as f64 + 0.13719) / (d as f64);
            Cpx::new(rad * ang.cos(), rad * ang.sin())
        })
        .collect();
    for _ in 0..300 {
        let mut mx = 0.0_f64;
        for i in 0..d {
            let z = roots[i];
            let mut val = Cpx::new(c[d], 0.0);
            for k in (0..d).rev() {
                val = val.mul(z).add(Cpx::new(c[k], 0.0));
            }
            let mut den = Cpx::new(1.0, 0.0);
            for j in 0..d {
                if i != j {
                    den = den.mul(z.sub(roots[j]));
                }
            }
            if den.abs() < 1e-24 { continue; }
            let dz = val.div(den);
            roots[i] = z.sub(dz);
            mx = mx.max(dz.abs());
        }
        if mx < 1e-12 { break; }
    }
    roots
}

fn spec_valid(sp: &Spec) -> bool {
    let m = sp.rs.len();
    if m < 2 || sp.ls.len() != m { return false; }
    let sum_l: i32 = sp.ls.iter().sum();
    if sum_l != sp.d1 + sp.d2 || sum_l <= 0 { return false; }
    if sp.d1 <= 0 || sp.d2 <= 0 || sp.s1 < 0 || sp.s2 < sp.s1 { return false; }
    // Ensure overlap of poles is non-empty and nested or overlapping cleanly
    let j_lo = sp.s2;
    let j_hi = (sp.s1 + sp.d1).min(sp.s2 + sp.d2);
    if j_lo >= j_hi { return false; }
    let e1 = sp.s1 + sp.d1 - sp.s2;
    let e2 = sp.s2 + sp.d2 - sp.s1;
    if e1 <= 0 || e2 <= 0 { return false; }
    for i in 0..m {
        if sp.ls[i] <= 0 { return false; }
        // require r_i - L_i >= 0 so all zeros are positive, and j + r_i - L_i >= 0 for all j >= s1
        if sp.rs[i] - sp.ls[i] < 0 { return false; }
    }
    // Check double-zero gap [T0, T1] is non-empty!
    let mut zmins: Vec<i32> = (0..m).map(|i| sp.rs[i] - sp.ls[i]).collect();
    zmins.sort();
    let mut rmaxs = sp.rs.clone();
    rmaxs.sort();
    let t0 = zmins[1];
    let t1 = rmaxs[m - 2];
    t1 > t0
}

#[derive(Clone, Debug)]
struct SaddleResult {
    sigma: f64,
    err_rate: f64,
    g: f64,
    t0: i32,
    m_max: i32,
}

fn eval_saddle(sp: &Spec) -> Option<SaddleResult> {
    if !spec_valid(sp) { return None; }
    let m = sp.rs.len();
    let mut p_l = vec![1.0_f64];
    for i in 0..m {
        p_l = poly_mul(&p_l, &[-(sp.rs[i] - sp.ls[i]) as f64, 1.0]);
    }
    p_l = poly_mul(&p_l, &[sp.s1 as f64, 1.0]);
    p_l = poly_mul(&p_l, &[sp.s2 as f64, 1.0]);

    let mut p_r = vec![1.0_f64];
    for i in 0..m {
        p_r = poly_mul(&p_r, &[-sp.rs[i] as f64, 1.0]);
    }
    p_r = poly_mul(&p_r, &[(sp.s1 + sp.d1) as f64, 1.0]);
    p_r = poly_mul(&p_r, &[(sp.s2 + sp.d2) as f64, 1.0]);

    let diff: Vec<f64> = p_l.iter().zip(p_r.iter()).map(|(a, b)| a - b).collect();
    let roots = durand_kerner(diff);

    let e1 = (sp.s1 + sp.d1 - sp.s2) as f64;
    let e2 = (sp.s2 + sp.d2 - sp.s1) as f64;
    let c_const = e1 * e1.ln() + e2 * e2.ln()
        - sp.ls.iter().map(|&l| { let lf = l as f64; lf * lf.ln() }).sum::<f64>();

    let j_lo = sp.s2 as f64;
    let j_hi = (sp.s1 + sp.d1).min(sp.s2 + sp.d2) as f64;

    let mut zmins: Vec<i32> = (0..m).map(|i| sp.rs[i] - sp.ls[i]).collect();
    zmins.sort();
    let mut rmaxs = sp.rs.clone();
    rmaxs.sort();
    let t0 = zmins[1];
    let t1 = rmaxs[m - 2];

    let mut pole_sigma = None;
    let mut best_err = f64::NEG_INFINITY;

    for z in &roots {
        if z.im.abs() < 1e-6 && z.re > -j_hi && z.re < -j_lo {
            let j = -z.re;
            let mut sig = c_const;
            for i in 0..m {
                let a = j + sp.rs[i] as f64;
                let b = j + (sp.rs[i] - sp.ls[i]) as f64;
                if a <= 0.0 || b < 0.0 { return None; }
                sig += a * a.ln() - if b > 1e-14 { b * b.ln() } else { 0.0 };
            }
            let d1a = j - sp.s1 as f64;
            let d1b = (sp.s1 + sp.d1) as f64 - j;
            let d2a = j - sp.s2 as f64;
            let d2b = (sp.s2 + sp.d2) as f64 - j;
            if d1a <= 0.0 || d1b <= 0.0 || d2a <= 0.0 || d2b <= 0.0 { return None; }
            sig -= d1a * d1a.ln() + d1b * d1b.ln() + d2a * d2a.ln() + d2b * d2b.ln();
            if pole_sigma.map_or(true, |cur| sig > cur) {
                pole_sigma = Some(sig);
            }
        } else if z.im > 1e-6 {
            let mut val = Cpx::new(c_const, 0.0);
            for i in 0..m {
                let t_a = Cpx::new(z.re - (sp.rs[i] - sp.ls[i]) as f64, z.im);
                let t_b = Cpx::new(z.re - sp.rs[i] as f64, z.im);
                val = val.add(t_a.xlogx()).sub(t_b.xlogx());
            }
            let p1a = Cpx::new(z.re + (sp.s1 + sp.d1) as f64, z.im);
            let p1b = Cpx::new(z.re + sp.s1 as f64, z.im);
            let p2a = Cpx::new(z.re + (sp.s2 + sp.d2) as f64, z.im);
            let p2b = Cpx::new(z.re + sp.s2 as f64, z.im);
            val = val.sub(p1a.xlogx()).add(p1b.xlogx()).sub(p2a.xlogx()).add(p2b.xlogx());
            let er = val.re - 2.0 * PI * z.im;
            if er > best_err { best_err = er; }
        } else if z.im.abs() < 1e-6 && z.re > (t0 as f64) - 0.5 && z.re < (t1 as f64) + 0.5 {
            // Real saddle inside double-zero gap!
            let mut val = c_const;
            for i in 0..m {
                let a = (z.re - (sp.rs[i] - sp.ls[i]) as f64).abs();
                let b = (z.re - sp.rs[i] as f64).abs();
                if a > 1e-14 { val += (z.re - (sp.rs[i] - sp.ls[i]) as f64) * a.ln(); }
                if b > 1e-14 { val -= (z.re - sp.rs[i] as f64) * b.ln(); }
            }
            for &(s, d) in &[(sp.s1, sp.d1), (sp.s2, sp.d2)] {
                let a = z.re + (s + d) as f64;
                let b = z.re + s as f64;
                if a <= 0.0 || b <= 0.0 { continue; }
                val -= a * a.ln() - b * b.ln();
            }
            if val > best_err { best_err = val; }
        }
    }

    let sigma = pole_sigma?;
    if !best_err.is_finite() { return None; }
    let g = sigma - best_err;
    if g <= 0.0 { return None; }
    let m_max = (sp.s1 + sp.d1).max(sp.s2 + sp.d2) + t0;
    Some(SaddleResult { sigma, err_rate: best_err, g, t0, m_max })
}

// Fast approximate carry integral via dense grid (for Phase 1 screening)
fn fast_carry_u_lt_1(sp: &Spec, t0: i32, m_max: i32) -> f64 {
    // Evaluate u in [1/m_max, 1.0] at N_grid Chebyshev/uniform points in 1/u (i.e. p/n in [1, m_max]!)
    // Since du / u^2 = -d(1/u), let w = p/n in [1, m_max]! Then integral is simply int_1^{m_max} v(1/w) dw!!
    // WHOA: in w = p/n, the measure is UNIFORM dw!!
    let e1 = (sp.s1 + sp.d1 - sp.s2) as f64;
    let e2 = (sp.s2 + sp.d2 - sp.s1) as f64;
    let max_den = sp.ls.iter().copied().max().unwrap().max(sp.s1 + sp.d1 - sp.s2).max(sp.s2 + sp.d2 - sp.s1) as f64;
    let s1 = sp.s1 as f64;
    let s2 = sp.s2 as f64;
    let p1 = (sp.s1 + sp.d1) as f64;
    let p2 = (sp.s2 + sp.d2) as f64;
    let t0f = t0 as f64;
    let m = sp.rs.len();

    let steps = 8 * (m_max as usize).max(16);
    let dw = (m_max as f64 - 1.0) / (steps as f64);
    let mut sum = 0.0;
    let mut tpts = Vec::with_capacity(48);

    for step in 0..steps {
        let w = 1.0 + (step as f64 + 0.5) * dw;
        let u = 1.0 / w;
        let j_lo = s1 * u;
        let j_hi = p1.max(p2) * u;
        tpts.clear();
        tpts.push(j_lo);
        tpts.push(j_hi);
        let a2 = s2 * u;
        let b2 = p2 * u;
        let b1 = p1 * u;
        if a2 > j_lo && a2 < j_hi { tpts.push(a2); }
        if b2 > j_lo && b2 < j_hi { tpts.push(b2); }
        if b1 > j_lo && b1 < j_hi { tpts.push(b1); }
        let th = 1.0 - t0f * u;
        if th > j_lo && th < j_hi { tpts.push(th); }

        // Add fractional jump points for t in [j_lo, j_hi]
        for i in 0..m {
            for &sh in &[sp.rs[i] as f64 * u, (sp.rs[i] - sp.ls[i]) as f64 * u] {
                let k0 = (j_lo + sh).floor() as i32;
                let k1 = (j_hi + sh).ceil() as i32;
                for k in k0..=k1 {
                    let c = k as f64 - sh;
                    if c > j_lo && c < j_hi { tpts.push(c); }
                }
            }
        }
        for &sh in &[-s1 * u, -s2 * u] {
            let k0 = (j_lo + sh).floor() as i32;
            let k1 = (j_hi + sh).ceil() as i32;
            for k in k0..=k1 {
                let c = k as f64 - sh;
                if c > j_lo && c < j_hi { tpts.push(c); }
            }
        }
        for &sh in &[p1 * u, p2 * u] {
            let k0 = (sh - j_hi).floor() as i32;
            let k1 = (sh - j_lo).ceil() as i32;
            for k in k0..=k1 {
                let c = sh - k as f64;
                if c > j_lo && c < j_hi { tpts.push(c); }
            }
        }
        tpts.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut mn = 999_i32;
        for idx in 0..tpts.len() - 1 {
            if tpts[idx + 1] - tpts[idx] < 1e-10 { continue; }
            let t = 0.5 * (tpts[idx] + tpts[idx + 1]);
            let mut vb = 0_i32;
            for i in 0..m {
                vb += (t + sp.rs[i] as f64 * u).floor() as i32
                    - (sp.ls[i] as f64 * u).floor() as i32
                    - (t + (sp.rs[i] - sp.ls[i]) as f64 * u).floor() as i32;
            }
            let in1 = t > s1 * u && t < p1 * u;
            let in2 = t > s2 * u && t < p2 * u;
            let val = if in1 && in2 {
                vb += (e1 * u).floor() as i32 - (t - s2 * u).floor() as i32 - (p1 * u - t).floor() as i32;
                vb += (e2 * u).floor() as i32 - (t - s1 * u).floor() as i32 - (p2 * u - t).floor() as i32;
                let h3 = if t + t0f * u >= 1.0 { 3 } else { 0 };
                let h2 = (if t + t0f * u >= 1.0 { 2 } else { 0 }) + (if max_den * u >= 1.0 { 1 } else { 0 });
                3 + vb - h3.max(h2)
            } else if in1 {
                if t < s2 * u {
                    vb += (e1 * u).floor() as i32 + (s2 * u - t).floor() as i32 - (p1 * u - t).floor() as i32;
                    vb += (e2 * u).floor() as i32 - (t - s1 * u).floor() as i32 - (p2 * u - t).floor() as i32;
                } else {
                    vb += (e1 * u).floor() as i32 - (t - s2 * u).floor() as i32 - (p1 * u - t).floor() as i32;
                    vb += (e2 * u).floor() as i32 + (t - p2 * u).floor() as i32 - (t - s1 * u).floor() as i32;
                }
                let h2 = if t + t0f * u >= 1.0 { 2 } else { 0 };
                3 + vb - h2
            } else {
                // in2 only (overlapping case where p2 > p1)
                if t < s1 * u {
                    vb += (e1 * u).floor() as i32 - (t - s2 * u).floor() as i32 - (p1 * u - t).floor() as i32;
                    vb += (e2 * u).floor() as i32 + (s1 * u - t).floor() as i32 - (p2 * u - t).floor() as i32;
                } else {
                    vb += (e1 * u).floor() as i32 + (t - p1 * u).floor() as i32 - (t - s2 * u).floor() as i32;
                    vb += (e2 * u).floor() as i32 - (t - s1 * u).floor() as i32 - (p2 * u - t).floor() as i32;
                }
                let h2 = if t + t0f * u >= 1.0 { 2 } else { 0 };
                3 + vb - h2
            };
            if val < mn { mn = val; }
        }
        sum += (mn as f64) * dw;
    }
    // Also add periodic u >= 1 estimate (~32 points on (0,1) weighted by trigamma(x) - 1/x^2)
    let p_steps = 64;
    let dx = 1.0 / (p_steps as f64);
    let mut ypts = Vec::with_capacity(24);
    for st in 0..p_steps {
        let x = (st as f64 + 0.5) * dx;
        ypts.clear();
        ypts.push(0.0);
        ypts.push(1.0);
        for i in 0..m {
            ypts.push((-(sp.rs[i] as f64) * x).rem_euclid(1.0));
            ypts.push((-((sp.rs[i] - sp.ls[i]) as f64) * x).rem_euclid(1.0));
        }
        for &c in &[s1, s2, p1, p2] {
            ypts.push((c * x).rem_euclid(1.0));
        }
        ypts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut mn = 999_i32;
        for idx in 0..ypts.len() - 1 {
            if ypts[idx + 1] - ypts[idx] < 1e-10 { continue; }
            let y = 0.5 * (ypts[idx] + ypts[idx + 1]);
            let mut v = 0_i32;
            for i in 0..m {
                v += (y + sp.rs[i] as f64 * x).floor() as i32
                    - (sp.ls[i] as f64 * x).floor() as i32
                    - (y + (sp.rs[i] - sp.ls[i]) as f64 * x).floor() as i32;
            }
            v += (e1 * x).floor() as i32 - (y - s2 * x).floor() as i32 - (p1 * x - y).floor() as i32;
            v += (e2 * x).floor() as i32 - (y - s1 * x).floor() as i32 - (p2 * x - y).floor() as i32;
            if v < mn { mn = v; }
        }
        if mn != 0 {
            // trigamma(x) - 1/x^2 = sum_{k=1}^inf 1/(x+k)^2
            let mut tg = 0.0;
            for k in 1..20 {
                let d = x + k as f64;
                tg += 1.0 / (d * d);
            }
            let z = x + 20.0;
            tg += 1.0 / z + 0.5 / (z * z) + 1.0 / (6.0 * z * z * z);
            sum += (mn as f64) * tg * dx;
        }
    }
    sum
}

// ============================================================================
// EXACT Z/p^4Z ORACLE (captures carry + ALL 1st-order & 2nd-order F_p shells!)
// ============================================================================

fn sieve_primes(limit: usize) -> Vec<u64> {
    let mut is_p = vec![true; limit + 1];
    is_p[0] = false;
    if limit >= 1 { is_p[1] = false; }
    let mut i = 2;
    while i * i <= limit {
        if is_p[i] {
            let mut j = i * i;
            while j <= limit { is_p[j] = false; j += i; }
        }
        i += 1;
    }
    (2..=limit).filter(|&k| is_p[k]).map(|k| k as u64).collect()
}

fn mod_inv_u128(a: u128, m: u128) -> u128 {
    let (mut t, mut newt) = (0_i128, 1_i128);
    let (mut r, mut newr) = (m as i128, (a % m) as i128);
    while newr != 0 {
        let q = r / newr;
        let tmp_t = t - q * newt;
        t = newt; newt = tmp_t;
        let tmp_r = r - q * newr;
        r = newr; newr = tmp_r;
    }
    if t < 0 { t += m as i128; }
    t as u128
}

// Exact p-adic saving per n for primes sqrt(N) < p <= M*n (plus exact periodic carry for u >= 1).
// Returns (total_saving_per_n, u_lt_1_saving_per_n, shell_bonus_over_carry_u_lt_1)
fn exact_padic_saving(sp: &Spec, n: i32, t0: i32, m_max: i32, primes: &[u64]) -> (f64, f64) {
    let n_u = n as usize;
    let e1 = ((sp.s1 + sp.d1 - sp.s2) * n) as usize;
    let e2 = ((sp.s2 + sp.d2 - sp.s1) * n) as usize;
    let max_r = *sp.rs.iter().max().unwrap();
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2);
    let max_arg = ((max_pole + max_r + t0) * n) as usize + 2;
    let p_max = (m_max * n) as u64;
    // Use p > sqrt(max_arg) AND p > n (focusing on u in (1/M, 1) where ALL F_p shells live,
    // wait: also include p in (sqrt(max_arg), n] so u >= 1 is ALSO measured exactly!)
    let p_min_sq = (max_arg as u64).max(30 * (n as u64));

    let jmin = (sp.s1.min(sp.s2) * n) as usize;
    let jmax = (max_pole * n) as usize;
    let s1n = (sp.s1 * n) as usize;
    let p1n = ((sp.s1 + sp.d1) * n) as usize;
    let s2n = (sp.s2 * n) as usize;
    let p2n = ((sp.s2 + sp.d2) * n) as usize;
    let t0n = (t0 * n) as usize;
    let sum_l_odd = (sp.ls.iter().sum::<i32>() * n) % 2 != 0;
    let m = sp.rs.len();

    let mut fv = vec![0_i32; max_arg + 1];
    let mut fu = vec![1_u128; max_arg + 1];
    let mut h1 = vec![0_u128; max_arg + 1];
    let mut h2 = vec![0_u128; max_arg + 1];
    let mut h3 = vec![0_u128; max_arg + 1];

    let mut tot_log = 0.0_f64;
    let mut u_lt_1_log = 0.0_f64;

    for &p in primes {
        if p > p_max { break; }
        if p * p <= p_min_sq { continue; }
        if p == 2 || p == 3 { continue; }
        // Use mod = p^4 (fits easily in u128 for p <= 50000: 50000^4 ~ 6e18 < 2^64, product < 2^128!)
        let p128 = p as u128;
        let p2 = p128 * p128;
        let p3 = p2 * p128;
        let mod_p = p3 * p128; // p^4; we can detect v_p up to 4!
        let inv2 = (mod_p + 1) / 2;

        fv[0] = 0; fu[0] = 1;
        h1[0] = 0; h2[0] = 0; h3[0] = 0;
        for k in 1..=max_arg {
            let ku = k as u128;
            if ku % p128 == 0 {
                let q = ku / p128; // coprime to p since p^2 > max_arg
                fv[k] = fv[k - 1] + 1;
                fu[k] = (fu[k - 1] * q) % mod_p;
                let invq = mod_inv_u128(q, mod_p);
                let invq2 = (invq * invq) % mod_p;
                let invq3 = (invq2 * invq) % mod_p;
                h1[k] = (h1[k - 1] + invq) % mod_p;
                h2[k] = (h2[k - 1] + invq2) % mod_p;
                h3[k] = (h3[k - 1] + invq3) % mod_p;
            } else {
                fv[k] = fv[k - 1];
                fu[k] = (fu[k - 1] * ku) % mod_p;
                let invk = mod_inv_u128(ku, mod_p);
                let invk2 = (invk * invk) % mod_p;
                let invk3 = (invk2 * invk) % mod_p;
                h1[k] = (h1[k - 1] + p128 * invk) % mod_p;
                h2[k] = (h2[k - 1] + p2 * invk2) % mod_p;
                h3[k] = (h3[k - 1] + p3 * invk3) % mod_p;
            }
        }

        // C'_n
        let mut cv = fv[e1] + fv[e2];
        let mut l_prod = 1_u128;
        for &l in &sp.ls {
            let ln = (l * n) as usize;
            cv -= fv[ln];
            l_prod = (l_prod * fu[ln]) % mod_p;
        }
        let cu = (((fu[e1] * fu[e2]) % mod_p) * mod_inv_u128(l_prod, mod_p)) % mod_p;
        let cu_signed = if sum_l_odd { (mod_p - cu) % mod_p } else { cu };

        let mut q_acc = 0_u128;
        let mut p_acc = 0_u128;

        for j in jmin..=jmax {
            let in1 = j >= s1n && j <= p1n;
            let in2 = j >= s2n && j <= p2n;
            if !in1 && !in2 { continue; }
            let mut nv = cv;
            let mut nu = cu_signed;
            let mut den_prod = 1_u128;
            let mut dlog_s1 = 0_u128;
            for i in 0..m {
                let a = j + (sp.rs[i] * n) as usize;
                let b = j + ((sp.rs[i] - sp.ls[i]) * n) as usize;
                nv += fv[a] - fv[b];
                nu = (nu * fu[a]) % mod_p;
                den_prod = (den_prod * fu[b]) % mod_p;
                let diff_h = (h1[a] + mod_p - h1[b]) % mod_p;
                dlog_s1 = (dlog_s1 + mod_p - diff_h) % mod_p;
            }
            let idx = j + t0n;
            if in1 && in2 {
                let a1 = j - s1n; let b1 = p1n - j;
                let a2 = j - s2n; let b2 = p2n - j;
                let dv = fv[a1] + fv[b1] + fv[a2] + fv[b2];
                let mut du = (((fu[a1] * fu[b1]) % mod_p) * ((fu[a2] * fu[b2]) % mod_p)) % mod_p;
                du = (du * den_prod) % mod_p;
                if (a1 + a2) % 2 == 1 { du = (mod_p - du) % mod_p; }
                let av = nv - dv;
                if av < 0 { return (0.0, 0.0); } // shouldn't happen
                if av < 4 {
                    let au = (nu * mod_inv_u128(du, mod_p)) % mod_p;
                    let dlog_den = (((h1[b1] + mod_p - h1[a1]) % mod_p) + ((h1[b2] + mod_p - h1[a2]) % mod_p)) % mod_p;
                    let diff_s1 = (dlog_s1 + mod_p - dlog_den) % mod_p;
                    let p_av = match av { 0 => 1, 1 => p128, 2 => p2, 3 => p3, _ => 0 };
                    if av == 0 {
                        q_acc = (q_acc + au * p3) % mod_p;
                    }
                    let t1 = (((au * h3[idx]) % mod_p) * p_av) % mod_p;
                    let t2 = (((((au * inv2) % mod_p) * diff_s1) % mod_p) * ((h2[idx] * p_av) % mod_p)) % mod_p;
                    p_acc = (p_acc + t1 + t2) % mod_p;
                }
            } else if in1 {
                let a1 = j - s1n; let b1 = p1n - j;
                let (dv, mut du, neg) = if j < s2n {
                    let m1 = p2n - j; let m0 = s2n - j - 1;
                    let dv = fv[a1] + fv[b1] + fv[m1] - fv[m0];
                    let num_extra = fu[m0];
                    nu = (nu * num_extra) % mod_p;
                    let du = (((fu[a1] * fu[b1]) % mod_p) * fu[m1]) % mod_p;
                    (dv, du, a1 % 2 == 1)
                } else {
                    let m1 = j - s2n; let m0 = j - p2n - 1;
                    let dv = fv[a1] + fv[b1] + fv[m1] - fv[m0];
                    nu = (nu * fu[m0]) % mod_p;
                    let du = (((fu[a1] * fu[b1]) % mod_p) * fu[m1]) % mod_p;
                    (dv, du, (a1 + (sp.d2 * n) as usize + 1) % 2 == 1)
                };
                du = (du * den_prod) % mod_p;
                if neg { du = (mod_p - du) % mod_p; }
                let bv = nv - dv;
                if bv + 1 < 0 { return (0.0, 0.0); }
                if bv + 1 < 4 {
                    let bu = (nu * mod_inv_u128(du, mod_p)) % mod_p;
                    let p_pow = match bv + 1 { 0 => 1, 1 => p128, 2 => p2, 3 => p3, _ => 0 };
                    let t2 = ((((bu * inv2) % mod_p) * h2[idx]) % mod_p * p_pow) % mod_p;
                    p_acc = (p_acc + t2) % mod_p;
                }
            } else {
                // in2 only
                let a2 = j - s2n; let b2 = p2n - j;
                let (dv, mut du, neg) = if j < s1n {
                    let m1 = p1n - j; let m0 = s1n - j - 1;
                    let dv = fv[a2] + fv[b2] + fv[m1] - fv[m0];
                    nu = (nu * fu[m0]) % mod_p;
                    let du = (((fu[a2] * fu[b2]) % mod_p) * fu[m1]) % mod_p;
                    (dv, du, a2 % 2 == 1)
                } else {
                    let m1 = j - s1n; let m0 = j - p1n - 1;
                    let dv = fv[a2] + fv[b2] + fv[m1] - fv[m0];
                    nu = (nu * fu[m0]) % mod_p;
                    let du = (((fu[a2] * fu[b2]) % mod_p) * fu[m1]) % mod_p;
                    (dv, du, (a2 + (sp.d1 * n) as usize + 1) % 2 == 1)
                };
                du = (du * den_prod) % mod_p;
                if neg { du = (mod_p - du) % mod_p; }
                let bv = nv - dv;
                if bv + 1 < 0 { return (0.0, 0.0); }
                if bv + 1 < 4 {
                    let bu = (nu * mod_inv_u128(du, mod_p)) % mod_p;
                    let p_pow = match bv + 1 { 0 => 1, 1 => p128, 2 => p2, 3 => p3, _ => 0 };
                    let t2 = ((((bu * inv2) % mod_p) * h2[idx]) % mod_p * p_pow) % mod_p;
                    p_acc = (p_acc + t2) % mod_p;
                }
            }
        }

        let val_mod = |mut x: u128| -> u32 {
            if x == 0 { return 4; }
            let mut v = 0;
            while x % p128 == 0 { x /= p128; v += 1; }
            v
        };
        let v = val_mod(q_acc).min(val_mod(p_acc));
        if v > 0 {
            let lg = (v as f64) * (p as f64).ln();
            tot_log += lg;
            if p > n as u64 {
                u_lt_1_log += lg;
            }
        }
    }
    let nf = n_u as f64;
    (tot_log / nf, u_lt_1_log / nf)
}

#[derive(Clone, Debug)]
struct Candidate {
    sp: Spec,
    g: f64,
    sigma: f64,
    m_max: i32,
    t0: i32,
    carry: f64,
    mu_car: f64,
    score_pot: f64,
}

fn eval_candidate(sp: Spec) -> Option<Candidate> {
    let sr = eval_saddle(&sp)?;
    let carry = fast_carry_u_lt_1(&sp, sr.t0, sr.m_max);
    let h_car = sr.sigma + (3 * sr.m_max) as f64 - carry;
    let gap_car = sr.g - h_car;
    let mu_car = if gap_car > 0.01 { sr.g / gap_car } else { 999.0 };
    // Shell-potential score: LK has carry_gap ~ 16.59, extra shells ~ 4.22 (~0.14*(D1+D2))
    let sum_l = (sp.d1 + sp.d2) as f64;
    let h_pot = h_car - 0.16 * sum_l;
    let gap_pot = sr.g - h_pot;
    let score_pot = if gap_pot > 0.01 { sr.g / gap_pot } else { 999.0 };
    Some(Candidate { sp, g: sr.g, sigma: sr.sigma, m_max: sr.m_max, t0: sr.t0, carry, mu_car, score_pot })
}

fn main() {
    let t_start = Instant::now();
    let n_threads = 16;
    let primes = Arc::new(sieve_primes(5000));

    // Benchmark LK first!
    let lk = Spec { rs: vec![10, 11, 12], ls: vec![10, 10, 10], s1: 1, d1: 16, s2: 2, d2: 14 };
    let lk_c = eval_candidate(lk.clone()).unwrap();
    let (lk_s30, lk_u30) = exact_padic_saving(&lk, 30, lk_c.t0, lk_c.m_max, &primes);
    let (lk_s55, lk_u55) = exact_padic_saving(&lk, 55, lk_c.t0, lk_c.m_max, &primes);
    let lk_h55 = lk_c.sigma + (3 * lk_c.m_max) as f64 - lk_s55;
    println!("=== BENCHMARK: Leibtag-Kaminer (10,11,12; 10,10,10; 1,16,2,14) ===");
    println!("  g = {:.6}, sigma = {:.6}, 3M = {}, carry = {:.6}, mu_car = {:.5}",
             lk_c.g, lk_c.sigma, 3 * lk_c.m_max, lk_c.carry, lk_c.mu_car);
    println!("  Exact p-adic n=30: sav={:.5} (u<1: {:.5}), n=55: sav={:.5} (u<1: {:.5}) -> emp mu(55) = {:.5} [paper asymp: 5.41966]",
             lk_s30, lk_u30, lk_s55, lk_u55, lk_c.g / (lk_c.g - lk_h55));

    // =========================================================================
    // PHASE 1: Systematic generation of 3-block (5F4) & 4-block (6F5) & 2-block
    // =========================================================================
    let mut specs: Vec<Spec> = Vec::new();

    // 1A. General 3-block 5F4 with equal L = L1=L2=L3 (LK's family generalized!)
    // L in 2..=14 -> D1+D2 = 3*L
    for l in 2..=14 {
        let tot = 3 * l;
        for d1 in (tot / 2)..=(tot - 2).min(2 * l + 4) {
            let d2 = tot - d1;
            if d2 < 2 { continue; }
            for s1 in 0..=4 {
                for s2 in s1..=(s1 + (d1 - d2) / 2 + 2).min(s1 + 5) {
                    // Shifts z_u = r_u - L in 0..=5
                    for z1 in 0..=3 {
                        for z2 in z1..=4 {
                            for z3 in z2..=6 {
                                let sp = Spec {
                                    rs: vec![z1 + l, z2 + l, z3 + l],
                                    ls: vec![l, l, l],
                                    s1, d1, s2, d2,
                                };
                                if spec_valid(&sp) { specs.push(sp); }
                            }
                        }
                    }
                }
            }
        }
    }

    // 1B. Asymmetric L1, L2, L3 3-block 5F4 (small/medium sizes + neighborhood of LK!)
    for l1 in 3..=12 {
        for l2 in l1..=13 {
            for l3 in l2..=14 {
                if l1 == l2 && l2 == l3 { continue; } // already in 1A
                if l3 - l1 > 4 { continue; }
                let tot = l1 + l2 + l3;
                if tot > 38 { continue; }
                for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 2) {
                    let d2 = tot - d1;
                    let s1 = 1;
                    let s2 = 1 + (d1 - d2) / 2;
                    for z1 in 0..=1 {
                        for z2 in z1..=2 {
                            for z3 in z2..=3 {
                                // Try permutations of (l1,l2,l3) assigned to (z1,z2,z3)
                                for &(a, b, c) in &[(l1,l2,l3), (l3,l2,l1), (l2,l1,l3)] {
                                    let sp = Spec {
                                        rs: vec![z1 + a, z2 + b, z3 + c],
                                        ls: vec![a, b, c],
                                        s1, d1, s2, d2,
                                    };
                                    if spec_valid(&sp) { specs.push(sp); }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 1C. 4-block 6F5 (equal L in 2..=10, and near-equal!)
    for l in 2..=10 {
        let tot = 4 * l;
        for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 2) {
            let d2 = tot - d1;
            for s1 in 0..=3 {
                for s2 in s1..=(s1 + (d1 - d2) / 2 + 2).min(s1 + 4) {
                    for z1 in 0..=2 {
                        for z2 in z1..=3 {
                            for z3 in z2..=4 {
                                for z4 in z3..=5 {
                                    let sp = Spec {
                                        rs: vec![z1 + l, z2 + l, z3 + l, z4 + l],
                                        ls: vec![l, l, l, l],
                                        s1, d1, s2, d2,
                                    };
                                    if spec_valid(&sp) { specs.push(sp); }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 1D. 5-block 7F6 (equal L in 2..=7)
    for l in 2..=7 {
        let tot = 5 * l;
        for d1 in (tot / 2)..=(tot / 2 + 3).min(tot - 2) {
            let d2 = tot - d1;
            let s1 = 1;
            let s2 = 1 + (d1 - d2) / 2;
            for z1 in 0..=1 {
                for z2 in z1..=2 {
                    for z3 in z2..=3 {
                        for z4 in z3..=4 {
                            for z5 in z4..=5 {
                                let sp = Spec {
                                    rs: vec![z1+l, z2+l, z3+l, z4+l, z5+l],
                                    ls: vec![l, l, l, l, l],
                                    s1, d1, s2, d2,
                                };
                                if spec_valid(&sp) { specs.push(sp); }
                            }
                        }
                    }
                }
            }
        }
    }

    // 1E. Scaled perturbations K * LK + delta for K = 1, 2
    for k in [1, 2] {
        for dr1 in -1..=1 {
        for dr2 in -1..=1 {
        for dr3 in -1..=1 {
        for dl in -1..=1 {
        for dd1 in -2..=2 {
        for ds1 in -1..=1 {
        for ds2 in -1..=1 {
            let l1 = 10 * k + dl;
            let l2 = 10 * k;
            let l3 = 10 * k - dl;
            let tot = l1 + l2 + l3;
            let d1 = 16 * k + dd1;
            let d2 = tot - d1;
            let sp = Spec {
                rs: vec![10*k + dr1, 11*k + dr2, 12*k + dr3],
                ls: vec![l1, l2, l3],
                s1: 1*k + ds1,
                d1,
                s2: 2*k + ds2,
                d2,
            };
            if spec_valid(&sp) { specs.push(sp); }
        }}}}}}}
    }

    // Deduplicate specs
    specs.sort_by(|a, b| (a.rs.len(), &a.rs, &a.ls, a.s1, a.d1, a.s2, a.d2)
        .cmp(&(b.rs.len(), &b.rs, &b.ls, b.s1, b.d1, b.s2, b.d2)));
    specs.dedup();
    println!("\n[Phase 1] Evaluating {} distinct m-block rational functions across {} threads...",
             specs.len(), n_threads);

    let specs_arc = Arc::new(specs);
    let results_mutex = Arc::new(Mutex::new(Vec::<Candidate>::new()));
    let chunk = (specs_arc.len() + n_threads - 1) / n_threads;

    let mut handles = Vec::new();
    for tid in 0..n_threads {
        let sp_ref = Arc::clone(&specs_arc);
        let res_ref = Arc::clone(&results_mutex);
        handles.push(thread::spawn(move || {
            let lo = tid * chunk;
            let hi = (lo + chunk).min(sp_ref.len());
            let mut local_best = Vec::new();
            for i in lo..hi {
                if let Some(c) = eval_candidate(sp_ref[i].clone()) {
                    if c.mu_car < 8.5 || c.score_pot < 6.0 {
                        local_best.push(c);
                    }
                }
            }
            let mut guard = res_ref.lock().unwrap();
            guard.extend(local_best);
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut cands = Arc::try_unwrap(results_mutex).unwrap().into_inner().unwrap();
    println!("[Phase 1] Done in {:.2}s. Retained {} promising candidates (mu_car<8.5 or pot<6.0).",
             t_start.elapsed().as_secs_f64(), cands.len());

    // Print best carry-only per block count m=3,4,5
    for m_cnt in [3, 4, 5] {
        let mut sub: Vec<&Candidate> = cands.iter().filter(|c| c.sp.rs.len() == m_cnt).collect();
        sub.sort_by(|a, b| a.mu_car.partial_cmp(&b.mu_car).unwrap());
        println!("\n--- Top 5 carry-only (m={}-block) ---", m_cnt);
        for c in sub.iter().take(5) {
            println!("  mu_car={:.5} pot={:.5} g={:.3} sig={:.3} 3M={} car={:.3} | r={:?} L={:?} s=({},{};{},{})",
                     c.mu_car, c.score_pot, c.g, c.sigma, 3*c.m_max, c.carry,
                     c.sp.rs, c.sp.ls, c.sp.s1, c.sp.d1, c.sp.s2, c.sp.d2);
        }
    }

    // Select diverse pool for Phase 2 exact p-adic shell evaluation!
    let mut pool: Vec<Candidate> = Vec::new();
    for m_cnt in [3, 4, 5] {
        let mut by_car: Vec<Candidate> = cands.iter().filter(|c| c.sp.rs.len() == m_cnt).cloned().collect();
        by_car.sort_by(|a, b| a.mu_car.partial_cmp(&b.mu_car).unwrap());
        pool.extend(by_car.into_iter().take(400));

        let mut by_pot: Vec<Candidate> = cands.iter().filter(|c| c.sp.rs.len() == m_cnt).cloned().collect();
        by_pot.sort_by(|a, b| a.score_pot.partial_cmp(&b.score_pot).unwrap());
        pool.extend(by_pot.into_iter().take(400));
    }
    // Always include LK
    pool.push(lk_c.clone());
    // Dedup pool
    pool.sort_by(|a, b| (a.sp.rs.len(), &a.sp.rs, &a.sp.ls, a.sp.s1, a.sp.d1, a.sp.s2, a.sp.d2)
        .cmp(&(b.sp.rs.len(), &b.sp.rs, &b.sp.ls, b.sp.s1, b.sp.d1, b.sp.s2, b.sp.d2)));
    pool.dedup_by(|a, b| a.sp == b.sp);
    cands.clear(); // free memory

    println!("\n[Phase 2] Running exact Z/p^4Z oracle (n=22) on {} candidates across {} threads...",
             pool.len(), n_threads);
    let t_p2 = Instant::now();
    let pool_arc = Arc::new(pool);
    let scored_mutex = Arc::new(Mutex::new(Vec::<(f64, f64, f64, Candidate)>::new()));
    let chunk2 = (pool_arc.len() + n_threads - 1) / n_threads;

    let mut handles2 = Vec::new();
    for tid in 0..n_threads {
        let p_ref = Arc::clone(&pool_arc);
        let pr_ref = Arc::clone(&primes);
        let sc_ref = Arc::clone(&scored_mutex);
        handles2.push(thread::spawn(move || {
            let lo = tid * chunk2;
            let hi = (lo + chunk2).min(p_ref.len());
            let mut local = Vec::new();
            for i in lo..hi {
                let c = &p_ref[i];
                // Keep p_max = m_max * n <= 2500
                let n_test = if c.m_max <= 25 { 24 } else if c.m_max <= 40 { 18 } else { 14 };
                let (sav, u1) = exact_padic_saving(&c.sp, n_test, c.t0, c.m_max, &pr_ref);
                if sav <= 0.0 { continue; }
                // Note: primes <= sqrt(max_arg) are skipped in exact_padic_saving;
                // roughly add periodic u>=1 carry minus skipped small primes, or compare directly via h_emp
                let h_emp = c.sigma + (3 * c.m_max) as f64 - sav;
                let mu_emp = if c.g > h_emp + 0.01 { c.g / (c.g - h_emp) } else { 999.0 };
                local.push((mu_emp, sav, u1, c.clone()));
            }
            sc_ref.lock().unwrap().extend(local);
        }));
    }
    for h in handles2 { h.join().unwrap(); }

    let mut scored = Arc::try_unwrap(scored_mutex).unwrap().into_inner().unwrap();
    scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("[Phase 2] Done in {:.2}s.", t_p2.elapsed().as_secs_f64());

    println!("\n=== TOP 25 OVERALL BY EXACT p-ADIC ORACLE (n~24) ===");
    for (rank, (mu_e, sav, u1, c)) in scored.iter().take(25).enumerate() {
        let bonus = u1 - c.carry;
        println!("#{:02} m={} mu_emp={:.5} (car={:.5}) | g={:.3} sig={:.3} 3M={} sav={:.3} (u<1={:.3}, sh_est={:+.3}) | r={:?} L={:?} s=({},{};{},{})",
                 rank+1, c.sp.rs.len(), mu_e, c.mu_car, c.g, c.sigma, 3*c.m_max, sav, u1, bonus,
                 c.sp.rs, c.sp.ls, c.sp.s1, c.sp.d1, c.sp.s2, c.sp.d2);
    }

    // =========================================================================
    // PHASE 3: High-n verification (n=55 and n=85!) on top 40 candidates + LK!
    // =========================================================================
    let mut top_specs: Vec<Candidate> = scored.iter().take(40).map(|x| x.3.clone()).collect();
    if !top_specs.iter().any(|c| c.sp == lk) { top_specs.push(lk_c.clone()); }
    println!("\n[Phase 3] High-n exact Z/p^4Z verification (n=60) on top {} candidates...", top_specs.len());
    let t_p3 = Instant::now();
    let top_arc = Arc::new(top_specs);
    let fin_mutex = Arc::new(Mutex::new(Vec::<(f64, f64, f64, Candidate)>::new()));
    let mut handles3 = Vec::new();
    let chunk3 = (top_arc.len() + n_threads - 1) / n_threads;
    for tid in 0..n_threads {
        let t_ref = Arc::clone(&top_arc);
        let pr_ref = Arc::clone(&primes);
        let f_ref = Arc::clone(&fin_mutex);
        handles3.push(thread::spawn(move || {
            let lo = tid * chunk3;
            let hi = (lo + chunk3).min(t_ref.len());
            let mut local = Vec::new();
            for i in lo..hi {
                let c = &t_ref[i];
                let n_hi = (4800 / (c.m_max.max(1) as usize)).min(65).max(30) as i32;
                let (sav, u1) = exact_padic_saving(&c.sp, n_hi, c.t0, c.m_max, &pr_ref);
                let h_emp = c.sigma + (3 * c.m_max) as f64 - sav;
                let mu_emp = if c.g > h_emp + 0.01 { c.g / (c.g - h_emp) } else { 999.0 };
                local.push((mu_emp, sav, n_hi as f64, c.clone()));
            }
            f_ref.lock().unwrap().extend(local);
        }));
    }
    for h in handles3 { h.join().unwrap(); }

    let mut fin = Arc::try_unwrap(fin_mutex).unwrap().into_inner().unwrap();
    fin.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("[Phase 3] Done in {:.2}s (total {:.2}s).", t_p3.elapsed().as_secs_f64(), t_start.elapsed().as_secs_f64());
    println!("\n=== TOP 20 VERIFIED AT HIGH n (n~60) ===");
    for (rank, (mu_e, sav, n_hi, c)) in fin.iter().take(20).enumerate() {
        let is_lk = if c.sp == lk { " <--- LEIBTAG-KAMINER" } else { "" };
        println!("#{:02} m={} n={:2} mu_emp={:.5} | g={:.4} sig={:.4} 3M={} sav={:.4} | r={:?} L={:?} s=({},{};{},{}) {}",
                 rank+1, c.sp.rs.len(), *n_hi as i32, mu_e, c.g, c.sigma, 3*c.m_max, sav,
                 c.sp.rs, c.sp.ls, c.sp.s1, c.sp.d1, c.sp.s2, c.sp.d2, is_lk);
    }
}
