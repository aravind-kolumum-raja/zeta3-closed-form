use std::collections::HashSet;
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

// Canonical signature of the RATIONAL FUNCTION F_n(t) up to integer shift t -> t + c*n!
// Notice F_n(t) is uniquely determined by:
// - multiset of root multiplicities on intervals (or simply sorted (r_u - L_u, r_u) shifted by s1,
//   EXCEPT wait: two different factorizations (r1,L1),(r2,L2) can give the SAME numerator polynomial
//   if their union of intervals has identical multiplicity function!).
// So the TRUE canonical key is the multiplicity vector of zeros (shifted by +s1) and poles (shifted by -s1)!
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct CanonKey {
    // Represented cleanly by sorted (start, end) of merged/exact multiplicity steps!
    // Even simpler: for integer intervals [a, b], compute net multiplicity on each unit interval [k, k+1)!
    zero_mults: Vec<(i32, i32)>, // (k, mult) for shifted zero intervals [z_u + s1, r_u + s1)
    pole_mults: Vec<(i32, i32)>, // (k, mult) for shifted pole intervals [s_i - s1, s_i + D_i - s1]
}

fn canonical_key(sp: &Spec) -> CanonKey {
    let shift = sp.s1;
    // Zeros of (t - r_u n)_{L_u n} are on (r_u - L_u, r_u]; under t' = t - s1*n?
    // Wait: if we shift t -> t + s1*n, then:
    //   poles (t + s1*n)_{D1*n+1} -> poles at [0, D1] and [s2-s1, s2+D2-s1]
    //   zeros (t - r_u*n)_{L_u*n} -> zeros at [r_u - L_u + s1, r_u + s1]!
    // Wait: if s1 == s2, we also ensure d1 >= d2 by canonicalizing poles via unit-interval mults!
    let mut z_min = i32::MAX;
    let mut z_max = i32::MIN;
    for i in 0..sp.rs.len() {
        z_min = z_min.min(sp.rs[i] - sp.ls[i] + shift);
        z_max = z_max.max(sp.rs[i] + shift);
    }
    let mut zm = Vec::new();
    for k in z_min..z_max {
        let mut c = 0;
        for i in 0..sp.rs.len() {
            if k >= sp.rs[i] - sp.ls[i] + shift && k < sp.rs[i] + shift { c += 1; }
        }
        if c > 0 { zm.push((k, c)); }
    }
    let p_max = sp.d1.max(sp.s2 - sp.s1 + sp.d2);
    let mut pm = Vec::new();
    for k in 0..p_max {
        let mut c = 0;
        if k < sp.d1 { c += 1; }
        if k >= sp.s2 - sp.s1 && k < sp.s2 - sp.s1 + sp.d2 { c += 1; }
        if c > 0 { pm.push((k, c)); }
    }
    CanonKey { zero_mults: zm, pole_mults: pm }
}

// Convert any Spec to its cleanest nested pole form (s1 <= s2, p2 <= p1) if possible, keeping same F_n(t)!
fn normalize_spec(mut sp: Spec) -> Spec {
    let p1 = sp.s1 + sp.d1;
    let p2 = sp.s2 + sp.d2;
    let lo = sp.s1.min(sp.s2);
    let hi = p1.max(p2);
    let in_lo = sp.s1.max(sp.s2);
    let in_hi = p1.min(p2);
    // Replace poles with nested [lo, hi] and [in_lo, in_hi] — EXACT same multiset of poles!
    sp.s1 = lo;
    sp.d1 = hi - lo;
    sp.s2 = in_lo;
    sp.d2 = in_hi - in_lo;
    // Shift so s1 = 1 (or keep s1 = 1 if all z_u >= 0 after shifting!)
    // Actually keep s1 >= 1 and z_u = r_u - L_u >= 0: shift by delta = max(0, 1 - s1) or min z_u!
    let min_z = (0..sp.rs.len()).map(|i| sp.rs[i] - sp.ls[i]).min().unwrap();
    // Shifting t -> t - delta*n increases s1, s2 by delta and decreases r_u by delta!
    // Put min_z = 0 (maximizes s1, keeps r_u minimal, clean!)
    if min_z > 0 {
        sp.s1 += min_z;
        sp.s2 += min_z;
        for r in &mut sp.rs { *r -= min_z; }
    }
    sp
}

#[derive(Clone, Copy, Debug)]
struct Cpx { re: f64, im: f64 }
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
        for (j, &y) in b.iter().enumerate() { r[i + j] += x * y; }
    }
    r
}

fn durand_kerner(mut coeffs: Vec<f64>) -> Vec<Cpx> {
    while coeffs.len() > 1 && coeffs.last().unwrap().abs() < 1e-10 { coeffs.pop(); }
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
            for k in (0..d).rev() { val = val.mul(z).add(Cpx::new(c[k], 0.0)); }
            let mut den = Cpx::new(1.0, 0.0);
            for j in 0..d { if i != j { den = den.mul(z.sub(roots[j])); } }
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
    let j_lo = sp.s2;
    let j_hi = (sp.s1 + sp.d1).min(sp.s2 + sp.d2);
    if j_lo >= j_hi { return false; }
    let e1 = sp.s1 + sp.d1 - sp.s2;
    let e2 = sp.s2 + sp.d2 - sp.s1;
    if e1 <= 0 || e2 <= 0 { return false; }
    for i in 0..m {
        if sp.ls[i] <= 0 || sp.rs[i] - sp.ls[i] < 0 { return false; }
    }
    let mut zmins: Vec<i32> = (0..m).map(|i| sp.rs[i] - sp.ls[i]).collect();
    zmins.sort();
    let mut rmaxs = sp.rs.clone();
    rmaxs.sort();
    rmaxs[m - 2] > zmins[1]
}

#[derive(Clone, Debug)]
struct Saddle { sigma: f64, g: f64, t0: i32, m_max: i32 }

fn eval_saddle(sp: &Spec) -> Option<Saddle> {
    if !spec_valid(sp) { return None; }
    let m = sp.rs.len();
    let mut p_l = vec![1.0_f64];
    for i in 0..m { p_l = poly_mul(&p_l, &[-(sp.rs[i] - sp.ls[i]) as f64, 1.0]); }
    p_l = poly_mul(&p_l, &[sp.s1 as f64, 1.0]);
    p_l = poly_mul(&p_l, &[sp.s2 as f64, 1.0]);

    let mut p_r = vec![1.0_f64];
    for i in 0..m { p_r = poly_mul(&p_r, &[-sp.rs[i] as f64, 1.0]); }
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
            let d1a = j - sp.s1 as f64; let d1b = (sp.s1 + sp.d1) as f64 - j;
            let d2a = j - sp.s2 as f64; let d2b = (sp.s2 + sp.d2) as f64 - j;
            if d1a <= 0.0 || d1b <= 0.0 || d2a <= 0.0 || d2b <= 0.0 { return None; }
            sig -= d1a * d1a.ln() + d1b * d1b.ln() + d2a * d2a.ln() + d2b * d2b.ln();
            if pole_sigma.map_or(true, |cur| sig > cur) { pole_sigma = Some(sig); }
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
    Some(Saddle { sigma, g, t0, m_max })
}

fn sieve_primes(limit: usize) -> Vec<u64> {
    let mut is_p = vec![true; limit + 1];
    is_p[0] = false;
    if limit >= 1 { is_p[1] = false; }
    let mut i = 2;
    while i * i <= limit {
        if is_p[i] { let mut j = i * i; while j <= limit { is_p[j] = false; j += i; } }
        i += 1;
    }
    (2..=limit).filter(|&k| is_p[k]).map(|k| k as u64).collect()
}

fn mod_inv_u128(a: u128, m: u128) -> u128 {
    let (mut t, mut newt) = (0_i128, 1_i128);
    let (mut r, mut newr) = (m as i128, (a % m) as i128);
    while newr != 0 {
        let q = r / newr;
        let tmp_t = t - q * newt; t = newt; newt = tmp_t;
        let tmp_r = r - q * newr; r = newr; newr = tmp_r;
    }
    if t < 0 { t += m as i128; }
    t as u128
}

fn vp_exact_single(sp: &Spec, n: usize, p: u64, t0: i32) -> i32 {
    let e1 = ((sp.s1 + sp.d1 - sp.s2) as usize) * n;
    let e2 = ((sp.s2 + sp.d2 - sp.s1) as usize) * n;
    let max_r = *sp.rs.iter().max().unwrap() as usize;
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2) as usize;
    let max_arg = (max_pole + max_r + t0 as usize) * n + 2;
    let p128 = p as u128;
    let p2 = p128 * p128;
    let p3 = p2 * p128;
    let p4 = p3 * p128;
    let mod_p = p4 * p128; // p^5
    let inv2 = (mod_p + 1) / 2;

    let mut fv = vec![0_i32; max_arg + 1];
    let mut fu = vec![1_u128; max_arg + 1];
    let mut h1 = vec![0_u128; max_arg + 1];
    let mut h2 = vec![0_u128; max_arg + 1];
    let mut h3 = vec![0_u128; max_arg + 1];

    for k in 1..=max_arg {
        let ku = k as u128;
        if ku % p128 == 0 {
            let q = ku / p128;
            fv[k] = fv[k - 1] + 1;
            fu[k] = (fu[k - 1] * q) % mod_p;
            let invq = mod_inv_u128(q, mod_p);
            let iq2 = (invq * invq) % mod_p;
            let iq3 = (iq2 * invq) % mod_p;
            h1[k] = (h1[k - 1] + invq) % mod_p;
            h2[k] = (h2[k - 1] + iq2) % mod_p;
            h3[k] = (h3[k - 1] + iq3) % mod_p;
        } else {
            fv[k] = fv[k - 1];
            fu[k] = (fu[k - 1] * ku) % mod_p;
            let invk = mod_inv_u128(ku, mod_p);
            let ik2 = (invk * invk) % mod_p;
            let ik3 = (ik2 * invk) % mod_p;
            h1[k] = (h1[k - 1] + p128 * invk) % mod_p;
            h2[k] = (h2[k - 1] + p2 * ik2) % mod_p;
            h3[k] = (h3[k - 1] + p3 * ik3) % mod_p;
        }
    }

    let mut cv = fv[e1] + fv[e2];
    let mut l_prod = 1_u128;
    for &l in &sp.ls {
        let ln = (l as usize) * n;
        cv -= fv[ln];
        l_prod = (l_prod * fu[ln]) % mod_p;
    }
    let cu = (((fu[e1] * fu[e2]) % mod_p) * mod_inv_u128(l_prod, mod_p)) % mod_p;
    let sum_l_odd = ((sp.ls.iter().sum::<i32>() as usize) * n) % 2 != 0;
    let cu_signed = if sum_l_odd { (mod_p - cu) % mod_p } else { cu };

    let s1n = (sp.s1 as usize) * n;
    let p1n = ((sp.s1 + sp.d1) as usize) * n;
    let s2n = (sp.s2 as usize) * n;
    let p2n = ((sp.s2 + sp.d2) as usize) * n;
    let t0n = (t0 as usize) * n;
    let jmin = s1n.min(s2n);
    let jmax = p1n.max(p2n);
    let m = sp.rs.len();

    let mut q_acc = 0_u128;
    let mut p_acc = 0_u128;
    let p_pows = [1_u128, p128, p2, p3, p4];

    for j in jmin..=jmax {
        let in1 = j >= s1n && j <= p1n;
        let in2 = j >= s2n && j <= p2n;
        if !in1 && !in2 { continue; }
        let mut nv = cv;
        let mut nu = cu_signed;
        let mut den_prod = 1_u128;
        let mut dlog_s1 = 0_u128;
        for i in 0..m {
            let a = j + (sp.rs[i] as usize) * n;
            let b = j + ((sp.rs[i] - sp.ls[i]) as usize) * n;
            nv += fv[a] - fv[b];
            nu = (nu * fu[a]) % mod_p;
            den_prod = (den_prod * fu[b]) % mod_p;
            dlog_s1 = (dlog_s1 + mod_p - ((h1[a] + mod_p - h1[b]) % mod_p)) % mod_p;
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
            if av < 0 { return -99; }
            if av < 5 {
                let au = (nu * mod_inv_u128(du, mod_p)) % mod_p;
                let dlog_den = (((h1[b1] + mod_p - h1[a1]) % mod_p) + ((h1[b2] + mod_p - h1[a2]) % mod_p)) % mod_p;
                let diff_s1 = (dlog_s1 + mod_p - dlog_den) % mod_p;
                if av + 3 < 5 { q_acc = (q_acc + au * p_pows[(av + 3) as usize]) % mod_p; }
                let pav = p_pows[av as usize];
                let t1 = (((au * h3[idx]) % mod_p) * pav) % mod_p;
                let t2 = (((((au * inv2) % mod_p) * diff_s1) % mod_p) * ((h2[idx] * pav) % mod_p)) % mod_p;
                p_acc = (p_acc + t1 + t2) % mod_p;
            }
        } else if in1 {
            let a1 = j - s1n; let b1 = p1n - j;
            let (dv, mut du, neg) = if j < s2n {
                let m1 = p2n - j; let m0 = s2n - j - 1;
                nu = (nu * fu[m0]) % mod_p;
                (fv[a1] + fv[b1] + fv[m1] - fv[m0], (((fu[a1] * fu[b1]) % mod_p) * fu[m1]) % mod_p, a1 % 2 == 1)
            } else {
                let m1 = j - s2n; let m0 = j - p2n - 1;
                nu = (nu * fu[m0]) % mod_p;
                (fv[a1] + fv[b1] + fv[m1] - fv[m0], (((fu[a1] * fu[b1]) % mod_p) * fu[m1]) % mod_p, (a1 + (sp.d2 as usize) * n + 1) % 2 == 1)
            };
            du = (du * den_prod) % mod_p;
            if neg { du = (mod_p - du) % mod_p; }
            let bv = nv - dv;
            if bv + 1 < 0 { return -99; }
            if bv + 1 < 5 {
                let bu = (nu * mod_inv_u128(du, mod_p)) % mod_p;
                let t2 = ((((bu * inv2) % mod_p) * h2[idx]) % mod_p * p_pows[(bv + 1) as usize]) % mod_p;
                p_acc = (p_acc + t2) % mod_p;
            }
        } else {
            let a2 = j - s2n; let b2 = p2n - j;
            let (dv, mut du, neg) = if j < s1n {
                let m1 = p1n - j; let m0 = s1n - j - 1;
                nu = (nu * fu[m0]) % mod_p;
                (fv[a2] + fv[b2] + fv[m1] - fv[m0], (((fu[a2] * fu[b2]) % mod_p) * fu[m1]) % mod_p, a2 % 2 == 1)
            } else {
                let m1 = j - s1n; let m0 = j - p1n - 1;
                nu = (nu * fu[m0]) % mod_p;
                (fv[a2] + fv[b2] + fv[m1] - fv[m0], (((fu[a2] * fu[b2]) % mod_p) * fu[m1]) % mod_p, (a2 + (sp.d1 as usize) * n + 1) % 2 == 1)
            };
            du = (du * den_prod) % mod_p;
            if neg { du = (mod_p - du) % mod_p; }
            let bv = nv - dv;
            if bv + 1 < 0 { return -99; }
            if bv + 1 < 5 {
                let bu = (nu * mod_inv_u128(du, mod_p)) % mod_p;
                let t2 = ((((bu * inv2) % mod_p) * h2[idx]) % mod_p * p_pows[(bv + 1) as usize]) % mod_p;
                p_acc = (p_acc + t2) % mod_p;
            }
        }
    }

    let val_mod = |mut x: u128| -> i32 {
        if x == 0 { return 5; }
        let mut v = 0;
        while x % p128 == 0 { x /= p128; v += 1; }
        v
    };
    val_mod(q_acc).min(val_mod(p_acc))
}

fn digamma_p1(z: f64) -> f64 {
    let mut s = 0.0;
    for k in 1..=20 { s -= 1.0 / (z + k as f64); }
    let w = z + 21.0;
    s + w.ln() - 0.5 / w - 1.0 / (12.0 * w * w) + 1.0 / (120.0 * w * w * w * w)
}

fn gcd_i32(mut a: i32, mut b: i32) -> i32 {
    while b != 0 { let t = a % b; a = b; b = t; }
    a
}

fn farey_seq(n: i32) -> Vec<(i32, i32)> {
    let mut pts = Vec::new();
    pts.push((0, 1));
    pts.push((1, 1));
    for d in 1..=n {
        for k in 1..d { if gcd_i32(k, d) == 1 { pts.push((k, d)); } }
    }
    pts.sort_by(|a, b| (a.0 as i64 * b.1 as i64).cmp(&(b.0 as i64 * a.1 as i64)));
    pts
}

// Ultra-paranoid exact Farey-partition certifier with `n_primes` independent primes per cell!
// Also checks disagreements across primes (if > 0, indicates a finer Farey split inside a cell!).
fn paranoid_farey_cert(sp: &Spec, sr: &Saddle, primes: &[u64], farey_boost: i32, n_primes: usize) -> Option<(f64, f64, f64, f64, usize)> {
    let max_r = *sp.rs.iter().max().unwrap();
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2);
    let max_arg_1 = (max_pole + max_r + sr.t0) as u64;
    let farey_n = (max_pole + max_r).max(sr.m_max) + farey_boost;
    let fpts = farey_seq(farey_n);
    let inv_m = 1.0 / (sr.m_max as f64);

    let mut d_u_lt_1 = 0.0;
    let mut d_u_ge_1 = 0.0;
    let mut disagreements = 0_usize;

    // 1. u in (1/M, 1)
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        if b <= inv_m + 1e-12 { continue; }
        let a_cl = a.max(inv_m);
        if b - a_cl < 1e-12 { continue; }
        let min_p = (((max_arg_1 as f64) * b).ceil() as u64 + 7).max(sr.m_max as u64 + 7);
        let mut vals = Vec::new();
        for &p in primes {
            if p < min_p || p > 2400 { continue; }
            let pf = p as f64;
            let margin = (b - a_cl) * 0.08;
            let lo = ((a_cl + margin) * pf).ceil() as i64;
            let hi = ((b - margin) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                // Alternate choosing near-lo, mid, near-hi across primes to detect any interior jump!
                let n = match vals.len() % 3 {
                    0 => ((lo + hi) / 2) as usize,
                    1 => lo as usize,
                    _ => hi as usize,
                };
                if (p as u64) * (p as u64) > max_arg_1 * (n as u64) + 5 {
                    let v = vp_exact_single(sp, n, p, sr.t0);
                    if v < 0 { return None; }
                    vals.push(v);
                    if vals.len() >= n_primes { break; }
                }
            }
        }
        if vals.is_empty() { return None; }
        let v_min = *vals.iter().min().unwrap();
        let v_max = *vals.iter().max().unwrap();
        if v_min != v_max { disagreements += 1; }
        if v_min > 0 { d_u_lt_1 += (v_min as f64) * (1.0 / a_cl - 1.0 / b); }
    }

    // 2. Periodic u >= 1
    let min_p_per = 2 * max_arg_1 + 11;
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        if b - a < 1e-12 { continue; }
        let mut vals = Vec::new();
        for &p in primes {
            if p < min_p_per || p > 2450 { continue; }
            let pf = p as f64;
            let margin = (b - a) * 0.08;
            let lo = ((a + margin) * pf).ceil() as i64;
            let hi = ((b - margin) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let k = ((lo + hi) / 2) as usize;
                let n = p as usize + k;
                if (p as u64) * (p as u64) > max_arg_1 * (n as u64) + 5 {
                    let v = vp_exact_single(sp, n, p, sr.t0);
                    if v < 0 { return None; }
                    vals.push(v);
                    if v == 0 || vals.len() >= 2 { break; }
                }
            }
        }
        if vals.is_empty() { return None; }
        let v_min = *vals.iter().min().unwrap();
        if v_min > 0 { d_u_ge_1 += (v_min as f64) * (digamma_p1(b) - digamma_p1(a)); }
    }

    let d_tot = d_u_lt_1 + d_u_ge_1;
    let h = sr.sigma + (3 * sr.m_max) as f64 - d_tot;
    let mu = if sr.g > h + 1e-6 { sr.g / (sr.g - h) } else { 999.0 };
    Some((d_u_lt_1, d_u_ge_1, d_tot, mu, disagreements))
}

// Fast unbiased asymptotic integrator in w = p/n in [1, M] using uniform w-bins!
fn fast_unbiased_mu(sp: &Spec, sr: &Saddle, primes: &[u64]) -> Option<(f64, f64)> {
    let max_r = *sp.rs.iter().max().unwrap();
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2);
    let max_arg_1 = (max_pole + max_r + sr.t0) as u64;
    // Uniform grid in w = p/n in [1, M]: 60 bins! Each bin width dw = (M-1)/60.
    // For each bin midpoint w_mid, u = 1/w_mid, pick smallest (n, p) with n/p close to u!
    let bins = 64;
    let dw = (sr.m_max as f64 - 1.0) / (bins as f64);
    let mut d_u_lt_1 = 0.0;
    for b_idx in 0..bins {
        let w = 1.0 + (b_idx as f64 + 0.5) * dw;
        let u = 1.0 / w;
        // Find prime p > max(max_arg_1 * u + 5, 41) where round(u*p)/p is very close to u
        let min_p = ((max_arg_1 as f64) * u).ceil() as u64 + 5;
        let min_p = min_p.max(sr.m_max as u64 + 3);
        let mut best: Option<(usize, u64, f64)> = None;
        for &p in primes {
            if p < min_p { continue; }
            if p > min_p + 180 { break; }
            let n = (u * (p as f64)).round() as usize;
            if n == 0 || n >= p as usize { continue; }
            let err = ((n as f64) / (p as f64) - u).abs();
            if best.map_or(true, |b| err < b.2) {
                best = Some((n, p, err));
                if err < 0.002 * u { break; }
            }
        }
        let (n, p, _) = best?;
        let v = vp_exact_single(sp, n, p, sr.t0);
        if v < 0 { return None; }
        d_u_lt_1 += (v as f64) * dw;
    }
    // Periodic u >= 1: 16 bins
    let p_per = *primes.iter().find(|&&q| q >= 2 * max_arg_1 + 11)?;
    let mut d_u_ge_1 = 0.0;
    let pbins = 16;
    for bi in 0..pbins {
        let a = (bi as f64) / (pbins as f64);
        let b = ((bi + 1) as f64) / (pbins as f64);
        let k = (((a + b) * 0.5) * (p_per as f64)).round() as usize;
        let k = k.max(1).min(p_per as usize - 1);
        let v = vp_exact_single(sp, p_per as usize + k, p_per, sr.t0);
        if v > 0 { d_u_ge_1 += (v as f64) * (digamma_p1(b) - digamma_p1(a)); }
    }
    let d_tot = d_u_lt_1 + d_u_ge_1;
    let h = sr.sigma + (3 * sr.m_max) as f64 - d_tot;
    let mu = if sr.g > h + 1e-6 { sr.g / (sr.g - h) } else { 999.0 };
    Some((mu, d_tot))
}

fn main() {
    let t0 = Instant::now();
    let primes = Arc::new(sieve_primes(2600));

    // Step 0: Ultra-paranoid 5-prime + Farey-boost verification of #03 vs LK!
    let lk = normalize_spec(Spec { rs: vec![10, 11, 12], ls: vec![10, 10, 10], s1: 1, d1: 16, s2: 2, d2: 14 });
    let c03 = normalize_spec(Spec { rs: vec![13, 15, 17], ls: vec![13, 13, 13], s1: 2, d1: 21, s2: 3, d2: 18 });
    println!("=== ULTRA-PARANOID 5-PRIME CHECK (with Farey order N+15) ===");
    for (name, sp) in [("LK", &lk), ("#03", &c03)] {
        let sr = eval_saddle(sp).unwrap();
        let (d1, d2, dt, mu, dis) = paranoid_farey_cert(sp, &sr, &primes, 15, 5).unwrap();
        let h = sr.sigma + (3 * sr.m_max) as f64 - dt;
        println!("  {}: mu={:.8} | g={:.8} sig={:.8} 3M={} D={:.8} ({:.6}+{:.6}) h={:.8} | cell_disagreements={}",
                 name, mu, sr.g, sr.sigma, 3*sr.m_max, dt, d1, d2, h, dis);
    }

    // Step 1: Build CANONICALLY DEDUPLICATED catalog of 3-block (5F4) & 4-block (6F5) rational functions!
    let mut seen: HashSet<CanonKey> = HashSet::new();
    let mut catalog: Vec<(Spec, Saddle)> = Vec::new();
    let mut try_add = |raw: Spec| {
        if !spec_valid(&raw) { return; }
        let sp = normalize_spec(raw);
        if !spec_valid(&sp) { return; }
        let key = canonical_key(&sp);
        if seen.contains(&key) { return; }
        if let Some(sr) = eval_saddle(&sp) {
            if sr.m_max <= 44 && sr.g - sr.sigma > 0.65 * (sr.m_max as f64) {
                seen.insert(key);
                catalog.push((sp, sr));
            }
        }
    };

    // 1A. Symmetric 3-block (L1=L2=L3=L) for L in 5..=20, general shifts (0, z2, z3), general poles (s1, d1; s2, d2)!
    for l in 5..=20 {
        let tot = 3 * l;
        for d1 in (tot / 2)..=(tot / 2 + 6).min(tot - 3) {
            let d2 = tot - d1;
            let ds = d1 - d2;
            // s1 in 1..=5, s2 around s1 + ds/2
            for s1 in 1..=5 {
                let s2_c = s1 + ds / 2;
                for s2 in (s2_c - 1).max(s1)..=(s2_c + 2).min(s1 + ds) {
                    for z2 in 1..=4 {
                        for z3 in (z2 + 1)..=(z2 + 5).min(8) {
                            try_add(Spec { rs: vec![l, l + z2, l + z3], ls: vec![l, l, l], s1, d1, s2, d2 });
                        }
                    }
                }
            }
        }
    }

    // 1B. Asymmetric 3-block (L1, L2, L3 not all equal!) around good aspect ratios (tot in 24..=45)!
    for l_mid in 8..=15 {
        for dl1 in -3..=3 {
            for dl3 in -3..=3 {
                if dl1 == 0 && dl3 == 0 { continue; }
                let l1 = l_mid + dl1;
                let l2 = l_mid;
                let l3 = l_mid + dl3;
                if l1 < 5 || l3 < 5 { continue; }
                let tot = l1 + l2 + l3;
                for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 3) {
                    let d2 = tot - d1;
                    let ds = d1 - d2;
                    for s1 in 1..=3 {
                        for s2 in (s1 + ds / 2).max(s1)..=(s1 + (ds + 1) / 2 + 1).min(s1 + ds) {
                            for z2 in 1..=3 {
                                for z3 in (z2 + 1)..=(z2 + 3) {
                                    try_add(Spec { rs: vec![l1, l2 + z2, l3 + z3], ls: vec![l1, l2, l3], s1, d1, s2, d2 });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 1C. 4-block 6F5 (equal L in 4..=12, general shifts & poles)!
    for l in 4..=12 {
        let tot = 4 * l;
        for d1 in (tot / 2)..=(tot / 2 + 5).min(tot - 3) {
            let d2 = tot - d1;
            let ds = d1 - d2;
            for s1 in 1..=4 {
                let s2 = s1 + ds / 2;
                for z2 in 1..=3 {
                    for z3 in (z2 + 1)..=(z2 + 3) {
                        for z4 in (z3 + 1)..=(z3 + 3).min(8) {
                            try_add(Spec { rs: vec![l, l+z2, l+z3, l+z4], ls: vec![l,l,l,l], s1, d1, s2, d2 });
                        }
                    }
                }
            }
        }
    }

    // 1D. 5-block 7F6 (L in 3..=8)
    for l in 3..=8 {
        let tot = 5 * l;
        for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 3) {
            let d2 = tot - d1;
            let ds = d1 - d2;
            for s1 in 1..=3 {
                let s2 = s1 + ds / 2;
                for dz in [vec![0,1,2,3,4], vec![0,1,2,3,5], vec![0,2,3,4,6], vec![0,2,4,6,8]] {
                    try_add(Spec {
                        rs: dz.iter().map(|&z| l + z).collect(),
                        ls: vec![l; 5], s1, d1, s2, d2
                    });
                }
            }
        }
    }

    println!("\n[Stage 1] Screening {} CANONICALLY DISTINCT rational functions across 16 threads...", catalog.len());
    let t_s1 = Instant::now();
    let cat_arc = Arc::new(catalog);
    let hits_m = Arc::new(Mutex::new(Vec::<(f64, f64, Spec, Saddle)>::new()));
    let n_threads = 16;
    let chunk = (cat_arc.len() + n_threads - 1) / n_threads;
    let mut handles = Vec::new();
    for tid in 0..n_threads {
        let c_ref = Arc::clone(&cat_arc);
        let p_ref = Arc::clone(&primes);
        let h_ref = Arc::clone(&hits_m);
        handles.push(thread::spawn(move || {
            let lo = tid * chunk;
            let hi = (lo + chunk).min(c_ref.len());
            let mut local = Vec::new();
            for i in lo..hi {
                let (ref sp, ref sr) = c_ref[i];
                if let Some((mu, dt)) = fast_unbiased_mu(sp, sr, &p_ref) {
                    if mu < 5.55 { local.push((mu, dt, sp.clone(), sr.clone())); }
                }
            }
            h_ref.lock().unwrap().extend(local);
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut stage1 = Arc::try_unwrap(hits_m).unwrap().into_inner().unwrap();
    stage1.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("[Stage 1] Done in {:.2}s. Found {} DISTINCT rational functions with fast-asymp mu < 5.55!",
             t_s1.elapsed().as_secs_f64(), stage1.len());

    // Stage 2: Exact Farey-partition certification (3 primes/cell) on TOP 50 DISTINCT winners!
    let mut top50: Vec<(Spec, Saddle)> = stage1.into_iter().take(50).map(|(_, _, sp, sr)| (sp, sr)).collect();
    if !top50.iter().any(|(sp, _)| canonical_key(sp) == canonical_key(&lk)) {
        top50.push((lk.clone(), eval_saddle(&lk).unwrap()));
    }
    println!("\n[Stage 2] EXACT Farey-partition certification (3 primes/cell, N+6) on top {} DISTINCT winners...", top50.len());
    let t_s2 = Instant::now();
    let top_arc = Arc::new(top50);
    let cert_m = Arc::new(Mutex::new(Vec::<(f64, f64, f64, f64, usize, Spec, Saddle)>::new()));
    let chunk2 = (top_arc.len() + n_threads - 1) / n_threads;
    let mut h2 = Vec::new();
    for tid in 0..n_threads {
        let t_ref = Arc::clone(&top_arc);
        let p_ref = Arc::clone(&primes);
        let c_ref = Arc::clone(&cert_m);
        h2.push(thread::spawn(move || {
            let lo = tid * chunk2;
            let hi = (lo + chunk2).min(t_ref.len());
            let mut local = Vec::new();
            for i in lo..hi {
                let (ref sp, ref sr) = t_ref[i];
                if let Some((d1, d2, dt, mu, dis)) = paranoid_farey_cert(sp, sr, &p_ref, 6, 3) {
                    local.push((mu, dt, d1, d2, dis, sp.clone(), sr.clone()));
                }
            }
            c_ref.lock().unwrap().extend(local);
        }));
    }
    for h in h2 { h.join().unwrap(); }

    let mut certs = Arc::try_unwrap(cert_m).unwrap().into_inner().unwrap();
    certs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("[Stage 2] Done in {:.2}s (total {:.2}s).", t_s2.elapsed().as_secs_f64(), t0.elapsed().as_secs_f64());

    println!("\n=== TOP 30 CANONICALLY DISTINCT RATIONAL FUNCTIONS BY EXACT FAREY ASYMPTOTIC mu(zeta(3)) ===");
    let lk_key = canonical_key(&lk);
    for (rank, (mu, dt, d1, d2, dis, sp, sr)) in certs.iter().take(30).enumerate() {
        let h = sr.sigma + (3 * sr.m_max) as f64 - dt;
        let tag = if canonical_key(sp) == lk_key { " <--- LEIBTAG-KAMINER (5.41965)" } else { "" };
        println!("#{:02} m={} mu={:.6} dis={} | g={:.4} sig={:.4} 3M={} D={:.4} ({:.3}+{:.3}) h={:.4} | r={:?} L={:?} s=({},{};{},{}) {}",
                 rank + 1, sp.rs.len(), mu, dis, sr.g, sr.sigma, 3*sr.m_max, dt, d1, d2, h,
                 sp.rs, sp.ls, sp.s1, sp.d1, sp.s2, sp.d2, tag);
    }
}
