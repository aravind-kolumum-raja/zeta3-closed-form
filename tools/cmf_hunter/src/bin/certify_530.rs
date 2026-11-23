use std::f64::consts::PI;
use std::time::Instant;

#[derive(Clone, Debug)]
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
    let lead = coeffs[d];
    let c: Vec<f64> = coeffs.iter().map(|&x| x / lead).collect();
    let rad = (1.0 + c[..d].iter().map(|x| x.abs()).fold(0.0_f64, f64::max)).min(60.0);
    let mut roots: Vec<Cpx> = (0..d).map(|k| {
        let ang = 2.0 * PI * (k as f64 + 0.13719) / (d as f64);
        Cpx::new(rad * ang.cos(), rad * ang.sin())
    }).collect();
    for _ in 0..400 {
        let mut mx = 0.0_f64;
        for i in 0..d {
            let z = roots[i];
            let mut val = Cpx::new(c[d], 0.0);
            for k in (0..d).rev() { val = val.mul(z).add(Cpx::new(c[k], 0.0)); }
            let mut den = Cpx::new(1.0, 0.0);
            for j in 0..d { if i != j { den = den.mul(z.sub(roots[j])); } }
            let dz = val.div(den);
            roots[i] = z.sub(dz);
            mx = mx.max(dz.abs());
        }
        if mx < 1e-14 { break; }
    }
    roots
}

fn eval_saddle(sp: &Spec) -> (f64, f64, i32, i32) {
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
    let mut zmins: Vec<i32> = (0..m).map(|i| sp.rs[i] - sp.ls[i]).collect();
    zmins.sort();
    let t0 = zmins[1];
    let m_max = (sp.s1 + sp.d1).max(sp.s2 + sp.d2) + t0;
    let j_lo = sp.s2 as f64;
    let j_hi = (sp.s1 + sp.d1).min(sp.s2 + sp.d2) as f64;
    let mut sigma = 0.0;
    let mut best_err = f64::NEG_INFINITY;
    for z in &roots {
        if z.im.abs() < 1e-7 && z.re > -j_hi && z.re < -j_lo {
            let j = -z.re;
            let mut sig = c_const;
            for i in 0..m {
                let a = j + sp.rs[i] as f64;
                let b = j + (sp.rs[i] - sp.ls[i]) as f64;
                sig += a * a.ln() - if b > 1e-14 { b * b.ln() } else { 0.0 };
            }
            let d1a = j - sp.s1 as f64; let d1b = (sp.s1 + sp.d1) as f64 - j;
            let d2a = j - sp.s2 as f64; let d2b = (sp.s2 + sp.d2) as f64 - j;
            sig -= d1a * d1a.ln() + d1b * d1b.ln() + d2a * d2a.ln() + d2b * d2b.ln();
            sigma = sig;
        } else if z.im > 1e-7 {
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
    (sigma, sigma - best_err, t0, m_max)
}

fn carry_u_lt_1_exact(sp: &Spec, t0: i32, u: f64) -> i32 {
    let e1 = (sp.s1 + sp.d1 - sp.s2) as f64;
    let e2 = (sp.s2 + sp.d2 - sp.s1) as f64;
    let max_den = sp.ls.iter().copied().max().unwrap()
        .max(sp.s1 + sp.d1 - sp.s2).max(sp.s2 + sp.d2 - sp.s1) as f64;
    let s1 = sp.s1 as f64; let s2 = sp.s2 as f64;
    let p1 = (sp.s1 + sp.d1) as f64; let p2 = (sp.s2 + sp.d2) as f64;
    let t0f = t0 as f64;
    let m = sp.rs.len();
    let j_lo = s1 * u; let j_hi = p1.max(p2) * u;
    let mut tpts = vec![j_lo, j_hi];
    for &x in &[s2 * u, p2 * u, p1 * u, 1.0 - t0f * u] {
        if x > j_lo && x < j_hi { tpts.push(x); }
    }
    for i in 0..m {
        for &sh in &[sp.rs[i] as f64 * u, (sp.rs[i] - sp.ls[i]) as f64 * u] {
            for k in ((j_lo + sh).floor() as i32)..=((j_hi + sh).ceil() as i32) {
                let c = k as f64 - sh;
                if c > j_lo && c < j_hi { tpts.push(c); }
            }
        }
    }
    for &sh in &[-s1 * u, -s2 * u] {
        for k in ((j_lo + sh).floor() as i32)..=((j_hi + sh).ceil() as i32) {
            let c = k as f64 - sh;
            if c > j_lo && c < j_hi { tpts.push(c); }
        }
    }
    for &sh in &[p1 * u, p2 * u] {
        for k in ((sh - j_hi).floor() as i32)..=((sh - j_lo).ceil() as i32) {
            let c = sh - k as f64;
            if c > j_lo && c < j_hi { tpts.push(c); }
        }
    }
    tpts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut mn = 999;
    for idx in 0..tpts.len() - 1 {
        if tpts[idx + 1] - tpts[idx] < 1e-11 { continue; }
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
        } else {
            if t < s2 * u {
                vb += (e1 * u).floor() as i32 + (s2 * u - t).floor() as i32 - (p1 * u - t).floor() as i32;
                vb += (e2 * u).floor() as i32 - (t - s1 * u).floor() as i32 - (p2 * u - t).floor() as i32;
            } else {
                vb += (e1 * u).floor() as i32 - (t - s2 * u).floor() as i32 - (p1 * u - t).floor() as i32;
                vb += (e2 * u).floor() as i32 + (t - p2 * u).floor() as i32 - (t - s1 * u).floor() as i32;
            }
            let h2 = if t + t0f * u >= 1.0 { 2 } else { 0 };
            3 + vb - h2
        };
        if val < mn { mn = val; }
    }
    mn
}

fn carry_per_exact(sp: &Spec, x: f64) -> i32 {
    let e1 = (sp.s1 + sp.d1 - sp.s2) as f64;
    let e2 = (sp.s2 + sp.d2 - sp.s1) as f64;
    let s1 = sp.s1 as f64; let s2 = sp.s2 as f64;
    let p1 = (sp.s1 + sp.d1) as f64; let p2 = (sp.s2 + sp.d2) as f64;
    let m = sp.rs.len();
    let mut ypts = vec![0.0, 1.0];
    for i in 0..m {
        ypts.push((-(sp.rs[i] as f64) * x).rem_euclid(1.0));
        ypts.push((-((sp.rs[i] - sp.ls[i]) as f64) * x).rem_euclid(1.0));
    }
    for &c in &[s1, s2, p1, p2] { ypts.push((c * x).rem_euclid(1.0)); }
    ypts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut mn = 999;
    for idx in 0..ypts.len() - 1 {
        if ypts[idx + 1] - ypts[idx] < 1e-11 { continue; }
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
    mn
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
            let iq = mod_inv_u128(q, mod_p);
            let iq2 = (iq * iq) % mod_p;
            let iq3 = (iq2 * iq) % mod_p;
            h1[k] = (h1[k - 1] + iq) % mod_p;
            h2[k] = (h2[k - 1] + iq2) % mod_p;
            h3[k] = (h3[k - 1] + iq3) % mod_p;
        } else {
            fv[k] = fv[k - 1];
            fu[k] = (fu[k - 1] * ku) % mod_p;
            let ik = mod_inv_u128(ku, mod_p);
            let ik2 = (ik * ik) % mod_p;
            let ik3 = (ik2 * ik) % mod_p;
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
    let mut q_acc = 0_u128;
    let mut p_acc = 0_u128;

    for j in s1n..=p1n {
        let in2 = j >= s2n && j <= p2n;
        let mut nv = cv;
        let mut nu = cu;
        let mut dp = 1_u128;
        let mut dlog_s1 = 0_u128;
        for i in 0..m {
            let a = j + (sp.rs[i] as usize) * n;
            let b = j + ((sp.rs[i] - sp.ls[i]) as usize) * n;
            nv += fv[a] - fv[b];
            nu = (nu * fu[a]) % mod_p;
            dp = (dp * fu[b]) % mod_p;
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

fn dissect(name: &str, sp: &Spec, farey_n: i32, primes: &[u64]) {
    let (sigma, g, t0, m_max) = eval_saddle(sp);
    let max_r = *sp.rs.iter().max().unwrap();
    let max_pole = (sp.s1 + sp.d1).max(sp.s2 + sp.d2);
    let max_arg_1 = (max_pole + max_r + t0) as u64;
    let fpts = farey_seq(farey_n);
    let inv_m = 1.0 / (m_max as f64);

    let mut d_car_u1 = 0.0;
    let mut d_sh1_u1 = 0.0; // excess == 1
    let mut d_sh2_u1 = 0.0; // excess >= 2 (2nd-order windows!)
    let mut shells: Vec<((i32,i32), (i32,i32), i32, i32, f64)> = Vec::new(); // (a, b, v_car, v_min, weight)

    for i in 0..fpts.len() - 1 {
        let fa = fpts[i]; let fb = fpts[i + 1];
        let a = fa.0 as f64 / fa.1 as f64;
        let b = fb.0 as f64 / fb.1 as f64;
        if b <= inv_m + 1e-12 { continue; }
        let a_cl = a.max(inv_m);
        if b - a_cl < 1e-12 { continue; }
        let xm = 0.5 * (a_cl + b);
        let vc = carry_u_lt_1_exact(sp, t0, xm);

        let min_p = (((max_arg_1 as f64) * b).ceil() as u64 + 7).max(m_max as u64 + 7);
        let mut vals = Vec::new();
        for &p in primes {
            if p < min_p || p > 2500 { continue; }
            let pf = p as f64;
            let margin = (b - a_cl) * 0.06;
            let lo = ((a_cl + margin) * pf).ceil() as i64;
            let hi = ((b - margin) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let n = match vals.len() % 3 {
                    0 => ((lo + hi) / 2) as usize,
                    1 => lo as usize,
                    _ => hi as usize,
                };
                if p * p > max_arg_1 * (n as u64) + 5 {
                    vals.push(vp_exact(sp, n, p, t0));
                    if vals.len() >= 10 { break; }
                }
            }
        }
        let vmin = *vals.iter().min().unwrap();
        let w = 1.0 / a_cl - 1.0 / b;
        d_car_u1 += (vc as f64) * w;
        let ex = vmin - vc;
        if ex == 1 { d_sh1_u1 += w; }
        else if ex >= 2 { d_sh1_u1 += w; d_sh2_u1 += ((ex - 1) as f64) * w; }
        if ex > 0 {
            shells.push((fa, fb, vc, vmin, (ex as f64) * w));
        }
    }

    // Periodic u >= 1
    let mut d_car_per = 0.0;
    let mut d_sh_per = 0.0;
    let min_p_per = 2 * max_arg_1 + 11;
    let mut per_shells = Vec::new();
    for i in 0..fpts.len() - 1 {
        let fa = fpts[i]; let fb = fpts[i + 1];
        let a = fa.0 as f64 / fa.1 as f64;
        let b = fb.0 as f64 / fb.1 as f64;
        let xm = 0.5 * (a + b);
        let vc = carry_per_exact(sp, xm);
        let mut vals = Vec::new();
        for &p in primes {
            if p < min_p_per || p > 2550 { continue; }
            let pf = p as f64;
            let margin = (b - a) * 0.06;
            let lo = ((a + margin) * pf).ceil() as i64;
            let hi = ((b - margin) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let k = ((lo + hi) / 2) as usize;
                let n = p as usize + k;
                if p * p > max_arg_1 * (n as u64) + 5 {
                    let v = vp_exact(sp, n, p, t0);
                    vals.push(v);
                    if v == vc || vals.len() >= 6 { break; }
                }
            }
        }
        let vmin = *vals.iter().min().unwrap();
        let w = digamma_p1(b) - digamma_p1(a);
        d_car_per += (vc as f64) * w;
        if vmin > vc {
            d_sh_per += ((vmin - vc) as f64) * w;
            per_shells.push((fa, fb, vc, vmin, ((vmin - vc) as f64) * w));
        }
    }

    let d_car = d_car_u1 + d_car_per;
    let d_1st = d_car + d_sh1_u1 + d_sh_per;
    let d_full = d_car + d_sh1_u1 + d_sh2_u1 + d_sh_per;
    let m3 = (3 * m_max) as f64;
    let mu_car = g / (g - (sigma + m3 - d_car));
    let mu_1st = g / (g - (sigma + m3 - d_1st));
    let mu_full = g / (g - (sigma + m3 - d_full));

    println!("\n======================================================================");
    println!("SURGICAL 10-PRIME DISSECTION: {} (Farey order N={})", name, farey_n);
    println!("  r={:?} L={:?} s=({},{};{},{}) | T0={} M={} 3M={}",
             sp.rs, sp.ls, sp.s1, sp.d1, sp.s2, sp.d2, t0, m_max, 3*m_max);
    println!("  g = {:.10} | sigma = {:.10}", g, sigma);
    println!("  Carry-only:    D_car  = {:.8} ({:.6} + {:.6}) => h = {:.6}, mu_car  = {:.6}",
             d_car, d_car_u1, d_car_per, sigma + m3 - d_car, mu_car);
    println!("  + 1st-order:   D_1st  = {:.8} (+sh1={:.6}, +per={:.6}) => h = {:.6}, mu_1st  = {:.6}",
             d_1st, d_sh1_u1, d_sh_per, sigma + m3 - d_1st, mu_1st);
    println!("  + 2nd-order:   D_full = {:.8} (+sh2={:.6}) => h = {:.6}, mu_full = {:.6}",
             d_full, d_sh2_u1, sigma + m3 - d_full, mu_full);

    // Merge adjacent shell intervals with identical excess (vmin - vc)!
    let mut merged: Vec<((i32,i32), (i32,i32), i32, f64)> = Vec::new();
    for (fa, fb, vc, vmin, wt) in shells {
        let ex = vmin - vc;
        if let Some(last) = merged.last_mut() {
            if last.1 == fa && last.2 == ex {
                last.1 = fb;
                last.3 += wt;
                continue;
            }
        }
        merged.push((fa, fb, ex, wt));
    }
    println!("  Merged u<1 F_p-shell intervals ({} total, verified across 10 primes each):", merged.len());
    for (fa, fb, ex, wt) in &merged {
        let ord = if *ex == 1 { "1st-order" } else { "2ND-ORDER!" };
        println!("    ({:2}/{:<2}, {:2}/{:<2})  excess=+{} [{}]  saving = {:.6}",
                 fa.0, fa.1, fb.0, fb.1, ex, ord, wt);
    }
    // Merge periodic shells too
    let mut pmerged: Vec<((i32,i32), (i32,i32), i32, f64)> = Vec::new();
    for (fa, fb, vc, vmin, wt) in per_shells {
        let ex = vmin - vc;
        if let Some(last) = pmerged.last_mut() {
            if last.1 == fa && last.2 == ex { last.1 = fb; last.3 += wt; continue; }
        }
        pmerged.push((fa, fb, ex, wt));
    }
    println!("  Merged periodic (u>=1) F_p-shell intervals ({} total):", pmerged.len());
    for (fa, fb, ex, wt) in &pmerged {
        println!("    ({:2}/{:<2}, {:2}/{:<2})  excess=+{}  saving = {:.6}",
                 fa.0, fa.1, fb.0, fb.1, ex, wt);
    }
}

fn main() {
    let t0 = Instant::now();
    let primes = sieve_primes(2600);
    let lk = Spec { rs: vec![10, 11, 12], ls: vec![10, 10, 10], s1: 1, d1: 16, s2: 2, d2: 14 };
    let c03 = Spec { rs: vec![13, 15, 17], ls: vec![13, 13, 13], s1: 2, d1: 21, s2: 3, d2: 18 };
    dissect("LEIBTAG-KAMINER (arXiv:2610.27799)", &lk, 29, &primes);
    dissect("NEW RECORD #01 (13,15,17; 13,13,13; 2,21; 3,18)", &c03, 40, &primes);
    dissect("NEW RECORD #01 CROSS-CHECK (finer Farey N=55)", &c03, 55, &primes);
    println!("\nDone in {:.2}s.", t0.elapsed().as_secs_f64());
}
