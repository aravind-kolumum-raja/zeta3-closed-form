use std::collections::HashSet;
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Spec { rs: Vec<i32>, ls: Vec<i32>, s1: i32, d1: i32, s2: i32, d2: i32 }

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
    fn xlogx(self) -> Self { if self.abs() < 1e-14 { Self::new(0.0, 0.0) } else { self.mul(self.ln()) } }
}

fn poly_mul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut r = vec![0.0; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() { for (j, &y) in b.iter().enumerate() { r[i + j] += x * y; } }
    r
}

fn durand_kerner(mut coeffs: Vec<f64>) -> Vec<Cpx> {
    while coeffs.len() > 1 && coeffs.last().unwrap().abs() < 1e-10 { coeffs.pop(); }
    let d = coeffs.len() - 1;
    if d == 0 { return vec![]; }
    let lead = coeffs[d];
    let c: Vec<f64> = coeffs.iter().map(|&x| x / lead).collect();
    let rad = (1.0 + c[..d].iter().map(|x| x.abs()).fold(0.0_f64, f64::max)).min(80.0);
    let mut roots: Vec<Cpx> = (0..d).map(|k| {
        let ang = 2.0 * PI * (k as f64 + 0.13719) / (d as f64);
        Cpx::new(rad * ang.cos(), rad * ang.sin())
    }).collect();
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

fn eval_saddle(sp: &Spec) -> Option<(f64, f64, i32, i32)> {
    let m = sp.rs.len();
    if sp.ls.iter().sum::<i32>() != sp.d1 + sp.d2 { return None; }
    if sp.d1 <= 0 || sp.d2 <= 0 || sp.s1 < 0 || sp.s2 < sp.s1 { return None; }
    if sp.s2 + sp.d2 > sp.s1 + sp.d1 { return None; } // enforce nested
    for i in 0..m { if sp.ls[i] <= 0 || sp.rs[i] - sp.ls[i] < 0 { return None; } }
    let mut zmins: Vec<i32> = (0..m).map(|i| sp.rs[i] - sp.ls[i]).collect();
    zmins.sort();
    let mut rmaxs = sp.rs.clone(); rmaxs.sort();
    let t0 = zmins[1]; let t1 = rmaxs[m - 2];
    if t1 <= t0 { return None; }
    let mut p_l = vec![1.0_f64];
    for i in 0..m { p_l = poly_mul(&p_l, &[-(sp.rs[i] - sp.ls[i]) as f64, 1.0]); }
    p_l = poly_mul(&p_l, &[sp.s1 as f64, 1.0]); p_l = poly_mul(&p_l, &[sp.s2 as f64, 1.0]);
    let mut p_r = vec![1.0_f64];
    for i in 0..m { p_r = poly_mul(&p_r, &[-sp.rs[i] as f64, 1.0]); }
    p_r = poly_mul(&p_r, &[(sp.s1 + sp.d1) as f64, 1.0]); p_r = poly_mul(&p_r, &[(sp.s2 + sp.d2) as f64, 1.0]);
    let diff: Vec<f64> = p_l.iter().zip(p_r.iter()).map(|(a, b)| a - b).collect();
    let roots = durand_kerner(diff);
    let e1 = (sp.s1 + sp.d1 - sp.s2) as f64; let e2 = (sp.s2 + sp.d2 - sp.s1) as f64;
    if e1 <= 0.0 || e2 <= 0.0 { return None; }
    let c_const = e1 * e1.ln() + e2 * e2.ln()
        - sp.ls.iter().map(|&l| { let lf = l as f64; lf * lf.ln() }).sum::<f64>();
    let mut sigma = None; let mut best_err = f64::NEG_INFINITY;
    for z in &roots {
        if z.im.abs() < 1e-6 && z.re > -((sp.s2 + sp.d2) as f64) && z.re < -(sp.s2 as f64) {
            let j = -z.re;
            let mut sig = c_const;
            for i in 0..m {
                let a = j + sp.rs[i] as f64; let b = j + (sp.rs[i] - sp.ls[i]) as f64;
                if a <= 0.0 || b < 0.0 { return None; }
                sig += a * a.ln() - if b > 1e-14 { b * b.ln() } else { 0.0 };
            }
            let d1a = j - sp.s1 as f64; let d1b = (sp.s1 + sp.d1) as f64 - j;
            let d2a = j - sp.s2 as f64; let d2b = (sp.s2 + sp.d2) as f64 - j;
            if d1a <= 0.0 || d1b <= 0.0 || d2a <= 0.0 || d2b <= 0.0 { return None; }
            sig -= d1a * d1a.ln() + d1b * d1b.ln() + d2a * d2a.ln() + d2b * d2b.ln();
            sigma = Some(sig);
        } else if z.im > 1e-6 {
            let mut val = Cpx::new(c_const, 0.0);
            for i in 0..m {
                val = val.add(Cpx::new(z.re - (sp.rs[i] - sp.ls[i]) as f64, z.im).xlogx())
                         .sub(Cpx::new(z.re - sp.rs[i] as f64, z.im).xlogx());
            }
            val = val.sub(Cpx::new(z.re + (sp.s1 + sp.d1) as f64, z.im).xlogx())
                     .add(Cpx::new(z.re + sp.s1 as f64, z.im).xlogx())
                     .sub(Cpx::new(z.re + (sp.s2 + sp.d2) as f64, z.im).xlogx())
                     .add(Cpx::new(z.re + sp.s2 as f64, z.im).xlogx());
            best_err = best_err.max(val.re - 2.0 * PI * z.im);
        }
    }
    let s = sigma?;
    if !best_err.is_finite() { return None; }
    let g = s - best_err;
    if g <= 0.0 { return None; }
    Some((s, g, t0, sp.s1 + sp.d1 + t0))
}

fn sieve_primes(limit: usize) -> Vec<u64> {
    let mut is_p = vec![true; limit + 1];
    is_p[0] = false; is_p[1] = false;
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
        let tt = t - q * newt; t = newt; newt = tt;
        let rr = r - q * newr; r = newr; newr = rr;
    }
    if t < 0 { t += m as i128; }
    t as u128
}

fn vp_exact(sp: &Spec, n: usize, p: u64, t0: i32) -> i32 {
    let e1 = ((sp.s1 + sp.d1 - sp.s2) as usize) * n;
    let e2 = ((sp.s2 + sp.d2 - sp.s1) as usize) * n;
    let max_r = *sp.rs.iter().max().unwrap() as usize;
    let max_pole = (sp.s1 + sp.d1) as usize;
    let max_arg = (max_pole + max_r + t0 as usize) * n + 2;
    let p128 = p as u128;
    let p2 = p128 * p128; let p3 = p2 * p128; let p4 = p3 * p128;
    let mod_p = p4 * p128; // p^5 (safe up to p <= 2500)
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
            let iq = mod_inv_u128(q, mod_p);
            let iq2 = (iq * iq) % mod_p; let iq3 = (iq2 * iq) % mod_p;
            h1[k] = (h1[k - 1] + iq) % mod_p;
            h2[k] = (h2[k - 1] + iq2) % mod_p;
            h3[k] = (h3[k - 1] + iq3) % mod_p;
        } else {
            fv[k] = fv[k - 1];
            fu[k] = (fu[k - 1] * ku) % mod_p;
            let ik = mod_inv_u128(ku, mod_p);
            let ik2 = (ik * ik) % mod_p; let ik3 = (ik2 * ik) % mod_p;
            h1[k] = (h1[k - 1] + p128 * ik) % mod_p;
            h2[k] = (h2[k - 1] + p2 * ik2) % mod_p;
            h3[k] = (h3[k - 1] + p3 * ik3) % mod_p;
        }
    }
    let mut cv = fv[e1] + fv[e2];
    let mut lp = 1_u128;
    for &l in &sp.ls { let ln = (l as usize) * n; cv -= fv[ln]; lp = (lp * fu[ln]) % mod_p; }
    let cu = (((fu[e1] * fu[e2]) % mod_p) * mod_inv_u128(lp, mod_p)) % mod_p;
    let cu = if ((sp.ls.iter().sum::<i32>() as usize) * n) % 2 != 0 { (mod_p - cu) % mod_p } else { cu };

    let s1n = (sp.s1 as usize) * n; let p1n = ((sp.s1 + sp.d1) as usize) * n;
    let s2n = (sp.s2 as usize) * n; let p2n = ((sp.s2 + sp.d2) as usize) * n;
    let t0n = (t0 as usize) * n;
    let m = sp.rs.len();
    let p_pows = [1_u128, p128, p2, p3, p4];
    let mut q_acc = 0_u128; let mut p_acc = 0_u128;

    for j in s1n..=p1n {
        let in2 = j >= s2n && j <= p2n;
        let mut nv = cv; let mut nu = cu; let mut dp = 1_u128; let mut dlog_s1 = 0_u128;
        for i in 0..m {
            let a = j + (sp.rs[i] as usize) * n;
            let b = j + ((sp.rs[i] - sp.ls[i]) as usize) * n;
            nv += fv[a] - fv[b];
            nu = (nu * fu[a]) % mod_p; dp = (dp * fu[b]) % mod_p;
            dlog_s1 = (dlog_s1 + mod_p - ((h1[a] + mod_p - h1[b]) % mod_p)) % mod_p;
        }
        let idx = j + t0n;
        let a1 = j - s1n; let b1 = p1n - j;
        if in2 {
            let a2 = j - s2n; let b2 = p2n - j;
            let dv = fv[a1] + fv[b1] + fv[a2] + fv[b2];
            let mut du = (((fu[a1] * fu[b1]) % mod_p) * ((fu[a2] * fu[b2]) % mod_p)) % mod_p;
            du = (du * dp) % mod_p;
            if (a1 + a2) % 2 == 1 { du = (mod_p - du) % mod_p; }
            let av = nv - dv;
            if av < 0 { return -99; }
            if av < 5 {
                let au = (nu * mod_inv_u128(du, mod_p)) % mod_p;
                let dld = (((h1[b1] + mod_p - h1[a1]) % mod_p) + ((h1[b2] + mod_p - h1[a2]) % mod_p)) % mod_p;
                let diff = (dlog_s1 + mod_p - dld) % mod_p;
                if av + 3 < 5 { q_acc = (q_acc + au * p_pows[(av + 3) as usize]) % mod_p; }
                let pav = p_pows[av as usize];
                let t1 = (((au * h3[idx]) % mod_p) * pav) % mod_p;
                let t2 = (((((au * inv2) % mod_p) * diff) % mod_p) * ((h2[idx] * pav) % mod_p)) % mod_p;
                p_acc = (p_acc + t1 + t2) % mod_p;
            }
        } else {
            let (dv, mut du, neg) = if j < s2n {
                let m1 = p2n - j; let m0 = s2n - j - 1;
                nu = (nu * fu[m0]) % mod_p;
                (fv[a1] + fv[b1] + fv[m1] - fv[m0], (((fu[a1] * fu[b1]) % mod_p) * fu[m1]) % mod_p, a1 % 2 == 1)
            } else {
                let m1 = j - s2n; let m0 = j - p2n - 1;
                nu = (nu * fu[m0]) % mod_p;
                (fv[a1] + fv[b1] + fv[m1] - fv[m0], (((fu[a1] * fu[b1]) % mod_p) * fu[m1]) % mod_p, (a1 + (sp.d2 as usize) * n + 1) % 2 == 1)
            };
            du = (du * dp) % mod_p;
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
    let vmod = |mut x: u128| -> i32 {
        if x == 0 { return 5; }
        let mut v = 0; while x % p128 == 0 { x /= p128; v += 1; } v
    };
    vmod(q_acc).min(vmod(p_acc))
}

fn digamma_p1(z: f64) -> f64 {
    let mut s = 0.0;
    for k in 1..=20 { s -= 1.0 / (z + k as f64); }
    let w = z + 21.0;
    s + w.ln() - 0.5 / w - 1.0 / (12.0 * w * w) + 1.0 / (120.0 * w * w * w * w)
}

fn gcd_i32(mut a: i32, mut b: i32) -> i32 { while b != 0 { let t = a % b; a = b; b = t; } a }

fn farey_seq(n: i32) -> Vec<(i32, i32)> {
    let mut pts = vec![(0, 1), (1, 1)];
    for d in 1..=n { for k in 1..d { if gcd_i32(k, d) == 1 { pts.push((k, d)); } } }
    pts.sort_by(|a, b| (a.0 as i64 * b.1 as i64).cmp(&(b.0 as i64 * a.1 as i64)));
    pts
}

fn fast_mu(sp: &Spec, primes: &[u64]) -> Option<(f64, f64, f64, i32)> {
    let (sigma, g, t0, m_max) = eval_saddle(sp)?;
    let max_r = *sp.rs.iter().max().unwrap() as u64;
    let max_pole = (sp.s1 + sp.d1) as u64;
    let max_arg_1 = max_pole + max_r + t0 as u64;
    if 2 * max_arg_1 + 20 > 2400 { return None; }
    let bins = 80;
    let dw = (m_max as f64 - 1.0) / (bins as f64);
    let mut d1 = 0.0;
    for bi in 0..bins {
        let w = 1.0 + (bi as f64 + 0.5) * dw;
        let u = 1.0 / w;
        let min_p = (((max_arg_1 as f64) * u).ceil() as u64 + 7).max(m_max as u64 + 5);
        let mut best: Option<(usize, u64, f64)> = None;
        for &p in primes {
            if p < min_p { continue; }
            if p > min_p + 200 || p > 2450 { break; }
            let n = (u * p as f64).round() as usize;
            if n == 0 || n >= p as usize { continue; }
            let err = ((n as f64) / (p as f64) - u).abs();
            if best.map_or(true, |b| err < b.2) {
                best = Some((n, p, err));
                if err < 0.0015 * u { break; }
            }
        }
        let (n, p, _) = best?;
        let v = vp_exact(sp, n, p, t0);
        if v < 0 { return None; }
        d1 += (v as f64) * dw;
    }
    let p_per = *primes.iter().find(|&&q| q >= 2 * max_arg_1 + 11 && q <= 2450)?;
    let mut d2 = 0.0;
    for bi in 0..16 {
        let a = bi as f64 / 16.0; let b = (bi + 1) as f64 / 16.0;
        let k = (((a + b) * 0.5) * p_per as f64).round() as usize;
        let k = k.max(1).min(p_per as usize - 1);
        let v = vp_exact(sp, p_per as usize + k, p_per, t0);
        if v > 0 { d2 += (v as f64) * (digamma_p1(b) - digamma_p1(a)); }
    }
    let h = sigma + (3 * m_max) as f64 - (d1 + d2);
    let mu = if g > h + 1e-6 { g / (g - h) } else { 999.0 };
    Some((mu, g, sigma, m_max))
}

fn exact_farey_mu(sp: &Spec, primes: &[u64]) -> Option<(f64, f64, f64, f64, i32)> {
    let (sigma, g, t0, m_max) = eval_saddle(sp)?;
    let max_r = *sp.rs.iter().max().unwrap();
    let max_pole = sp.s1 + sp.d1;
    let max_arg_1 = (max_pole + max_r + t0) as u64;
    let farey_n = (max_pole + max_r).max(m_max);
    let fpts = farey_seq(farey_n);
    let inv_m = 1.0 / (m_max as f64);
    let mut d1 = 0.0;
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        if b <= inv_m + 1e-12 { continue; }
        let a_cl = a.max(inv_m);
        if b - a_cl < 1e-12 { continue; }
        let min_p = (((max_arg_1 as f64) * b).ceil() as u64 + 7).max(m_max as u64 + 7);
        let mut vmin = 99; let mut cnt = 0;
        for &p in primes {
            if p < min_p || p > 2480 { continue; }
            let pf = p as f64; let mrg = (b - a_cl) * 0.08;
            let lo = ((a_cl + mrg) * pf).ceil() as i64;
            let hi = ((b - mrg) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let n = ((lo + hi) / 2) as usize;
                if p * p > max_arg_1 * (n as u64) + 5 {
                    let v = vp_exact(sp, n, p, t0);
                    if v < 0 { return None; }
                    if v < vmin { vmin = v; }
                    cnt += 1; if cnt >= 3 { break; }
                }
            }
        }
        if cnt == 0 { return None; }
        if vmin > 0 { d1 += (vmin as f64) * (1.0 / a_cl - 1.0 / b); }
    }
    let min_p_per = 2 * max_arg_1 + 11;
    let mut d2 = 0.0;
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        let mut vmin = 99; let mut cnt = 0;
        for &p in primes {
            if p < min_p_per || p > 2480 { continue; }
            let pf = p as f64; let mrg = (b - a) * 0.08;
            let lo = ((a + mrg) * pf).ceil() as i64;
            let hi = ((b - mrg) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let k = ((lo + hi) / 2) as usize;
                let n = p as usize + k;
                if p * p > max_arg_1 * (n as u64) + 5 {
                    let v = vp_exact(sp, n, p, t0);
                    if v < 0 { return None; }
                    if v < vmin { vmin = v; }
                    cnt += 1; if vmin == 0 || cnt >= 2 { break; }
                }
            }
        }
        if cnt == 0 { return None; }
        if vmin > 0 { d2 += (vmin as f64) * (digamma_p1(b) - digamma_p1(a)); }
    }
    let dt = d1 + d2;
    let h = sigma + (3 * m_max) as f64 - dt;
    let mu = g / (g - h);
    Some((mu, g, sigma, dt, m_max))
}

fn main() {
    let t0 = Instant::now();
    let primes = Arc::new(sieve_primes(2500));
    let mut cands: Vec<Spec> = Vec::new();
    let mut seen: HashSet<Spec> = HashSet::new();
    let mut add = |sp: Spec| {
        if eval_saddle(&sp).is_some() && seen.insert(sp.clone()) { cands.push(sp); }
    };

    // 1. Symmetric arithmetic-progression families L in 12..=30 with step/L near 0.10..0.22!
    for l in 12..=30 {
        let tot = 3 * l;
        for step in 1..=6 {
            let ratio = step as f64 / l as f64;
            if ratio < 0.08 || ratio > 0.24 { continue; }
            // Near-equal step perturbations (z2, z3) around (step, 2*step)
            for dz2 in -1..=1 {
                for dz3 in -1..=1 {
                    let z2 = step + dz2; let z3 = 2 * step + dz3;
                    if z2 <= 0 || z3 <= z2 { continue; }
                    for d1 in (tot / 2 + 1)..=(tot / 2 + step + 2) {
                        let d2 = tot - d1;
                        if d2 < 2 { continue; }
                        let ds = d1 - d2;
                        for s1 in 1..=(step + 2) {
                            for s2 in (s1 + 1)..=(s1 + ds - 1).max(s1 + 1) {
                                if s2 + d2 <= s1 + d1 {
                                    add(Spec { rs: vec![l, l + z2, l + z3], ls: vec![l, l, l], s1, d1, s2, d2 });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 2. Direct fine perturbations around 2 * #01 = r=[26,30,34], L=[26,26,26], s=(4,42;6,36)!
    for dl1 in -2..=2 {
    for dl2 in -2..=2 {
    for dr2 in -1..=1 {
    for dr3 in -1..=1 {
    for dd in -2..=2 {
    for ds1 in -1..=1 {
    for ds2 in -1..=1 {
        let l1 = 26 + dl1; let l2 = 26 + dl2; let l3 = 26 - dl1 - dl2;
        let d1 = 42 + dd; let d2 = 78 - d1;
        let s1 = 4 + ds1; let s2 = 6 + ds2;
        add(Spec { rs: vec![l1, l2 + 4 + dr2, l3 + 8 + dr3], ls: vec![l1, l2, l3], s1, d1, s2, d2 });
    }}}}}}}

    println!("[Probe L=12..30 & 2x#01] Screening {} candidates across 16 threads...", cands.len());
    let c_arc = Arc::new(cands);
    let hits_m = Arc::new(Mutex::new(Vec::<(f64, Spec)>::new()));
    let n_threads = 16;
    let chunk = (c_arc.len() + n_threads - 1) / n_threads;
    let mut handles = Vec::new();
    for tid in 0..n_threads {
        let cr = Arc::clone(&c_arc); let pr = Arc::clone(&primes); let hr = Arc::clone(&hits_m);
        handles.push(thread::spawn(move || {
            let lo = tid * chunk; let hi = (lo + chunk).min(cr.len());
            let mut local = Vec::new();
            for i in lo..hi {
                if let Some((mu, _, _, _)) = fast_mu(&cr[i], &pr) {
                    if mu < 5.42 { local.push((mu, cr[i].clone())); }
                }
            }
            hr.lock().unwrap().extend(local);
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut hits = Arc::try_unwrap(hits_m).unwrap().into_inner().unwrap();
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("Done screening in {:.2}s. Found {} with fast_mu < 5.42. Running EXACT Farey cert on top 20...",
             t0.elapsed().as_secs_f64(), hits.len());

    let top20: Vec<Spec> = hits.into_iter().take(20).map(|x| x.1).collect();
    let t_arc = Arc::new(top20);
    let fin_m = Arc::new(Mutex::new(Vec::<(f64, f64, f64, f64, i32, Spec)>::new()));
    let mut h2 = Vec::new();
    let ch2 = (t_arc.len() + n_threads - 1) / n_threads;
    for tid in 0..n_threads {
        let tr = Arc::clone(&t_arc); let pr = Arc::clone(&primes); let fr = Arc::clone(&fin_m);
        h2.push(thread::spawn(move || {
            let lo = tid * ch2; let hi = (lo + ch2).min(tr.len());
            let mut local = Vec::new();
            for i in lo..hi {
                if let Some((mu, g, sig, dt, mm)) = exact_farey_mu(&tr[i], &pr) {
                    local.push((mu, g, sig, dt, mm, tr[i].clone()));
                }
            }
            fr.lock().unwrap().extend(local);
        }));
    }
    for h in h2 { h.join().unwrap(); }
    let mut fin = Arc::try_unwrap(fin_m).unwrap().into_inner().unwrap();
    fin.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("\n=== TOP 15 HIGH-L / 2x#01 CANDIDATES BY EXACT FAREY ASYMPTOTIC mu(zeta(3)) ===");
    for (rk, (mu, g, sig, dt, mm, sp)) in fin.iter().take(15).enumerate() {
        let h = sig + (3 * mm) as f64 - dt;
        println!("#{:02} mu={:.6} | g={:.4} sig={:.4} 3M={} D={:.4} h={:.4} | r={:?} L={:?} s=({},{};{},{})",
                 rk+1, mu, g, sig, 3*mm, dt, h, sp.rs, sp.ls, sp.s1, sp.d1, sp.s2, sp.d2);
    }
    println!("Total elapsed: {:.2}s.", t0.elapsed().as_secs_f64());
}
