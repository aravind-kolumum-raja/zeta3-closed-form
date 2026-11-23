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

// Compute exact v_p(p^3 * C'_n * (Q_n, P_{3,n})) for a single (n, p) with p^2 > max_arg, p <= M*n.
// Works mod p^5 (detects v up to 5!).
fn vp_exact_single(sp: &Spec, n: usize, p: u64, t0: i32) -> i32 {
    let e1 = ((sp.s1 + sp.d1 - sp.s2) as usize) * n;
    let e2 = ((sp.s2 + sp.d2 - sp.s1) as usize) * n;
    let max_r = *sp.rs.iter().max().unwrap() as usize;
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2) as usize;
    let max_arg = (max_pole + max_r + t0 as usize) * n + 2;
    let p128 = p as u128;
    // mod = p^5 (safe in u128 for p <= 2500: 2500^5 = 9.7e16, sq = 9.5e33 < 3.4e38)
    let p2 = p128 * p128;
    let p3 = p2 * p128;
    let p4 = p3 * p128;
    let mod_p = p4 * p128;
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
                if av + 3 < 5 {
                    q_acc = (q_acc + au * p_pows[(av + 3) as usize]) % mod_p;
                }
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

// Build Farey sequence of order N in [0, 1]
fn farey_seq(n: i32) -> Vec<(i32, i32)> {
    let mut pts = Vec::new();
    pts.push((0, 1));
    pts.push((1, 1));
    for d in 1..=n {
        for k in 1..d {
            if gcd_i32(k, d) == 1 { pts.push((k, d)); }
        }
    }
    pts.sort_by(|a, b| (a.0 as i64 * b.1 as i64).cmp(&(b.0 as i64 * a.1 as i64)));
    pts
}

// Find smallest (n, p) with p prime, p > min_p, a < n/p < b!
fn find_sample_n_p(a: f64, b: f64, min_p: u64, primes: &[u64], shift_1: bool) -> Option<(usize, u64)> {
    for &p in primes {
        if p < min_p { continue; }
        let pf = p as f64;
        // Want n/p strictly inside (a, b) with tiny margin
        let margin = (b - a) * 0.05;
        let lo = ((a + margin) * pf).ceil() as i64;
        let hi = ((b - margin) * pf).floor() as i64;
        if lo <= hi && lo >= 1 {
            let k = ((lo + hi) / 2) as usize;
            let n = if shift_1 { p as usize + k } else { k };
            return Some((n, p));
        }
    }
    None
}

// Exact Farey-partition asymptotic D'_n (BOTH u in (1/M, 1) AND periodic u >= 1)!
// Returns (D_u_lt_1, D_u_ge_1, D_tot, mu_asymp, shell_cells_count)
fn exact_farey_asymptotic(sp: &Spec, sr: &Saddle, primes: &[u64], verbose: bool) -> Option<(f64, f64, f64, f64, usize)> {
    let max_r = *sp.rs.iter().max().unwrap();
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2);
    let max_arg_1 = (max_pole + max_r + sr.t0) as u64;
    // Farey order covers all breakpoints!
    let farey_n = (max_pole + max_r).max(sr.m_max);
    let fpts = farey_seq(farey_n);
    let inv_m = 1.0 / (sr.m_max as f64);

    let mut d_u_lt_1 = 0.0;
    let mut d_u_ge_1 = 0.0;
    let mut active_cells = 0_usize;

    // 1. u in (1/M, 1)
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        if b <= inv_m + 1e-12 { continue; }
        let a_cl = a.max(inv_m);
        if b - a_cl < 1e-12 { continue; }
        // Need p^2 > max_arg_1 * n -> since n < b * p, p > max_arg_1 * b + 3!
        let min_p = ((max_arg_1 as f64) * b).ceil() as u64 + 5;
        let min_p = min_p.max(sr.m_max as u64 + 5);
        // Sample TWO distinct primes to guard against accidental congruence!
        let mut samples = Vec::new();
        for &p in primes {
            if p < min_p || p > 2200 { continue; }
            let pf = p as f64;
            let margin = (b - a_cl) * 0.08;
            let lo = ((a_cl + margin) * pf).ceil() as i64;
            let hi = ((b - margin) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let n = ((lo + hi) / 2) as usize;
                if (p as u64) * (p as u64) > max_arg_1 * (n as u64) + 5 {
                    samples.push((n, p));
                    if samples.len() >= 2 { break; }
                }
            }
        }
        if samples.is_empty() { return None; }
        let mut v_min = 99;
        for &(n, p) in &samples {
            let v = vp_exact_single(sp, n, p, sr.t0);
            if v < 0 { return None; }
            if v < v_min { v_min = v; }
        }
        if v_min > 0 {
            let w = 1.0 / a_cl - 1.0 / b;
            d_u_lt_1 += (v_min as f64) * w;
            active_cells += 1;
            if verbose && v_min >= 1 {
                // Only print interesting cells where v_min differs from trivial
                // (we'll print summary below)
            }
        }
    }

    // 2. Periodic u >= 1 (sample at n = p + k, so u = 1 + k/p in (1, 2); need p^2 > 2 * max_arg_1 * p => p > 2 * max_arg_1 + 5!)
    let min_p_per = 2 * max_arg_1 + 11;
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        if b - a < 1e-12 { continue; }
        // Fast precheck: is periodic carry or valuation possibly > 0? Just sample 1 prime (or 2 if > 0)!
        let mut v_min = 99;
        let mut cnt = 0;
        for &p in primes {
            if p < min_p_per || p > 2400 { continue; }
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
                    if v < v_min { v_min = v; }
                    cnt += 1;
                    if v_min == 0 || cnt >= 2 { break; }
                }
            }
        }
        if cnt == 0 { return None; }
        if v_min > 0 {
            d_u_ge_1 += (v_min as f64) * (digamma_p1(b) - digamma_p1(a));
        }
    }

    let d_tot = d_u_lt_1 + d_u_ge_1;
    let h = sr.sigma + (3 * sr.m_max) as f64 - d_tot;
    let mu = if sr.g > h + 1e-6 { sr.g / (sr.g - h) } else { 999.0 };
    if verbose {
        println!("  [EXACT FAREY N={}] 3M={} | D(u<1)={:.8} + D(u>=1)={:.8} = {:.8} | 3M-D={:.8} | h={:.8} | g={:.8} => mu={:.8}",
                 farey_n, 3*sr.m_max, d_u_lt_1, d_u_ge_1, d_tot, (3*sr.m_max) as f64 - d_tot, h, sr.g, mu);
    }
    Some((d_u_lt_1, d_u_ge_1, d_tot, mu, active_cells))
}

// Fast approximate Farey-midpoint asymptotic D'_n using a SINGLE prime p ~ 211 (for searching!)
fn fast_asymptotic_mu(sp: &Spec, sr: &Saddle, primes: &[u64]) -> Option<(f64, f64)> {
    // Integrate in w = p/n in [1, M] using ~40 exact samples across (1/M, 1) + 20 periodic samples!
    let max_r = *sp.rs.iter().max().unwrap();
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2);
    let max_arg_1 = (max_pole + max_r + sr.t0) as u64;
    // Choose a prime p >= max(83, max_arg_1 + 10)
    let p_target = (max_arg_1 + 10).max(97);
    let p = *primes.iter().find(|&&q| q >= p_target)?;
    let pf = p as f64;
    let inv_m = 1.0 / (sr.m_max as f64);
    let m = sr.m_max;
    // Evaluate EVERY n in [floor(p/M)+1 .. p-1] plus fractional edge cell!
    // Wait: near u = 1/M, put a single test point at exact midpoint of (1/M, n_start/p) using a larger prime!
    let n_start = ((p as i32) / m) + 1;
    let mut d_u_lt_1 = 0.0;
    // Edge interval (1/M, (n_start - 0.5)/p) if nonempty
    let first_b = ((n_start as f64) - 0.5) / pf;
    if first_b > inv_m + 1e-6 {
        if let Some((ne, pe)) = find_sample_n_p(inv_m, first_b, max_arg_1 + 5, primes, false) {
            let v = vp_exact_single(sp, ne, pe, sr.t0);
            if v > 0 { d_u_lt_1 += (v as f64) * (1.0 / inv_m - 1.0 / first_b); }
        }
    }
    for n in (n_start as usize)..(p as usize) {
        let a = (((n as f64) - 0.5) / pf).max(inv_m);
        let b = (((n as f64) + 0.5) / pf).min(1.0);
        if b <= a { continue; }
        let v = vp_exact_single(sp, n, p, sr.t0);
        if v < 0 { return None; }
        if v > 0 { d_u_lt_1 += (v as f64) * (1.0 / a - 1.0 / b); }
    }
    // Quick periodic u>=1 estimate using 24 evenly spaced points with p_per > 2*max_arg_1 + 10
    let p_per = *primes.iter().find(|&&q| q >= 2 * max_arg_1 + 11)?;
    let mut d_u_ge_1 = 0.0;
    let bins = 24;
    for b_idx in 0..bins {
        let a = (b_idx as f64) / (bins as f64);
        let b = ((b_idx + 1) as f64) / (bins as f64);
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
    let primes = Arc::new(sieve_primes(3000));

    // 1. CALIBRATE ON LEIBTAG-KAMINER & #03!
    let lk = Spec { rs: vec![10, 11, 12], ls: vec![10, 10, 10], s1: 1, d1: 16, s2: 2, d2: 14 };
    println!("=== CALIBRATION 1: Leibtag-Kaminer (10,11,12; 10,10,10; 1,16; 2,14) ===");
    let lk_sr = eval_saddle(&lk).unwrap();
    exact_farey_asymptotic(&lk, &lk_sr, &primes, true);

    let c03 = Spec { rs: vec![13, 15, 17], ls: vec![13, 13, 13], s1: 2, d1: 21, s2: 3, d2: 18 };
    println!("\n=== TEST #03: (13,15,17; 13,13,13; 2,21; 3,18) ===");
    let c03_sr = eval_saddle(&c03).unwrap();
    exact_farey_asymptotic(&c03, &c03_sr, &primes, true);

    // 2. SYSTEMATIC EXACT-ASYMPTOTIC SEARCH AROUND LK + GENERAL 3-BLOCK & 4-BLOCK!
    // Generate rich pool of symmetric & asymmetric 3-block (5F4) and 4-block (6F5) candidates
    let mut seeds: Vec<Spec> = Vec::new();
    // 2A. All symmetric 3-block with L in 4..=16 and small shift patterns!
    for l in 4..=16 {
        let tot = 3 * l;
        for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 2) {
            let d2 = tot - d1;
            let s1 = 1;
            let s2_base: i32 = 1 + (d1 - d2) / 2;
            for s2 in (s2_base.saturating_sub(1).max(1))..=(s2_base + 1) {
                for z1 in 0..=1 {
                    for dz1 in 1..=3 {
                        for dz2 in 1..=3 {
                            let z2 = z1 + dz1;
                            let z3 = z2 + dz2;
                            if z3 > 6 { continue; }
                            let sp = Spec { rs: vec![z1+l, z2+l, z3+l], ls: vec![l, l, l], s1, d1, s2, d2 };
                            if spec_valid(&sp) { seeds.push(sp); }
                        }
                    }
                }
            }
        }
    }
    // 2B. Asymmetric perturbations around LK (r=[10,11,12], L=[10,10,10], s=(1,16;2,14)) AND around L=11,12,13!
    for base_l in [9, 10, 11, 12, 13] {
        for dl1 in -2..=2 {
        for dl2 in -2..=2 {
        for dr1 in -1..=1 {
        for dr2 in -1..=1 {
        for dr3 in -1..=1 {
        for dd in -2..=2 {
            let l1 = base_l + dl1;
            let l2 = base_l + dl2;
            let l3 = base_l - dl1 - dl2;
            if l1 < 4 || l2 < 4 || l3 < 4 { continue; }
            let tot = l1 + l2 + l3;
            let d1 = (tot / 2) + 1 + dd;
            let d2 = tot - d1;
            let s1 = 1;
            let s2 = 2;
            let sp = Spec { rs: vec![l1 + dr1, l2 + 1 + dr2, l3 + 2 + dr3], ls: vec![l1, l2, l3], s1, d1, s2, d2 };
            if spec_valid(&sp) { seeds.push(sp); }
        }}}}}}
    }
    // 2C. 4-block 6F5 symmetric with L in 4..=9
    for l in 4..=9 {
        let tot = 4 * l;
        for d1 in (tot / 2)..=(tot / 2 + 3) {
            let d2 = tot - d1;
            let s1 = 1;
            let s2 = 1 + (d1 - d2) / 2;
            for z2 in 1..=2 {
                for z3 in z2..=3 {
                    for z4 in z3..=4 {
                        let sp = Spec { rs: vec![l, z2+l, z3+l, z4+l], ls: vec![l,l,l,l], s1, d1, s2, d2 };
                        if spec_valid(&sp) { seeds.push(sp); }
                    }
                }
            }
        }
    }
    // 2D. Also 2-block 4F3 (RV01 family!) for comparison
    for l1 in 6..=16 {
        for l2 in l1..=18 {
            let tot = l1 + l2;
            for d1 in (tot/2)..=(tot/2 + 3) {
                let d2 = tot - d1;
                for z2 in 0..=3 {
                    let sp = Spec { rs: vec![l1, l2+z2], ls: vec![l1, l2], s1: 1, d1, s2: 1 + (d1-d2)/2, d2 };
                    if spec_valid(&sp) { seeds.push(sp); }
                }
            }
        }
    }

    // Deduplicate & pre-filter by fast saddle check (g > 20 and g / (g - sigma) reasonable)
    seeds.sort_by(|a, b| (a.rs.len(), &a.rs, &a.ls, a.s1, a.d1, a.s2, a.d2)
        .cmp(&(b.rs.len(), &b.rs, &b.ls, b.s1, b.d1, b.s2, b.d2)));
    seeds.dedup();
    let filtered: Vec<(Spec, Saddle)> = seeds.into_iter().filter_map(|sp| {
        let sr = eval_saddle(&sp)?;
        // Quick heuristic filter: LK has sigma/g = 0.580, 3M/g = 0.479
        if sr.g - sr.sigma < (0.65 * sr.m_max as f64) { return None; }
        Some((sp, sr))
    }).collect();
    println!("\n[Stage A] Fast asymptotic evaluation on {} filtered candidates across 16 threads...", filtered.len());

    let filt_arc = Arc::new(filtered);
    let results_m = Arc::new(Mutex::new(Vec::<(f64, f64, Spec, Saddle)>::new()));
    let n_threads = 16;
    let chunk = (filt_arc.len() + n_threads - 1) / n_threads;
    let mut handles = Vec::new();
    for tid in 0..n_threads {
        let f_ref = Arc::clone(&filt_arc);
        let p_ref = Arc::clone(&primes);
        let r_ref = Arc::clone(&results_m);
        handles.push(thread::spawn(move || {
            let lo = tid * chunk;
            let hi = (lo + chunk).min(f_ref.len());
            let mut local = Vec::new();
            for i in lo..hi {
                let (ref sp, ref sr) = f_ref[i];
                if let Some((mu, dtot)) = fast_asymptotic_mu(sp, sr, &p_ref) {
                    if mu < 5.75 {
                        local.push((mu, dtot, sp.clone(), sr.clone()));
                    }
                }
            }
            r_ref.lock().unwrap().extend(local);
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut stage_a = Arc::try_unwrap(results_m).unwrap().into_inner().unwrap();
    stage_a.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("[Stage A] Done in {:.2}s. Found {} candidates with fast-asymp mu < 5.75.",
             t0.elapsed().as_secs_f64(), stage_a.len());

    // Stage B: Exact Farey-partition asymptotic evaluation on top 35 unique candidates!
    let mut top_b: Vec<(Spec, Saddle)> = stage_a.into_iter().take(40).map(|(_, _, sp, sr)| (sp, sr)).collect();
    if !top_b.iter().any(|(sp, _)| sp == &lk) { top_b.push((lk.clone(), lk_sr.clone())); }
    println!("\n[Stage B] Running EXACT Farey-partition asymptotic certifier on top {} candidates...", top_b.len());

    let top_arc = Arc::new(top_b);
    let exact_m = Arc::new(Mutex::new(Vec::<(f64, f64, f64, f64, Spec, Saddle)>::new()));
    let chunk_b = (top_arc.len() + n_threads - 1) / n_threads;
    let mut hb = Vec::new();
    for tid in 0..n_threads {
        let t_ref = Arc::clone(&top_arc);
        let p_ref = Arc::clone(&primes);
        let e_ref = Arc::clone(&exact_m);
        hb.push(thread::spawn(move || {
            let lo = tid * chunk_b;
            let hi = (lo + chunk_b).min(t_ref.len());
            let mut local = Vec::new();
            for i in lo..hi {
                let (ref sp, ref sr) = t_ref[i];
                if let Some((d1, d2, dtot, mu, _)) = exact_farey_asymptotic(sp, sr, &p_ref, false) {
                    local.push((mu, dtot, d1, d2, sp.clone(), sr.clone()));
                }
            }
            e_ref.lock().unwrap().extend(local);
        }));
    }
    for h in hb { h.join().unwrap(); }

    let mut exact_res = Arc::try_unwrap(exact_m).unwrap().into_inner().unwrap();
    exact_res.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    println!("\n=== TOP 25 BY EXACT FAREY-PARTITION ASYMPTOTIC mu(zeta(3)) ===");
    for (rank, (mu, dtot, d1, d2, sp, sr)) in exact_res.iter().take(25).enumerate() {
        let h = sr.sigma + (3 * sr.m_max) as f64 - dtot;
        let tag = if sp == &lk { " <--- LEIBTAG-KAMINER (paper: 5.41966)" } else { "" };
        println!("#{:02} m={} mu={:.6} | g={:.5} sig={:.5} 3M={} D={:.5} ({:.4}+{:.4}) h={:.5} | r={:?} L={:?} s=({},{};{},{}) {}",
                 rank + 1, sp.rs.len(), mu, sr.g, sr.sigma, 3*sr.m_max, dtot, d1, d2, h,
                 sp.rs, sp.ls, sp.s1, sp.d1, sp.s2, sp.d2, tag);
    }
    println!("\nTotal elapsed: {:.2}s.", t0.elapsed().as_secs_f64());
}
