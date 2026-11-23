use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

// ============================================================================
// Double-Double (DD) arithmetic: ~31 decimal digits (106 bits mantissa)
// Zero heap allocations, pure stack registers, < 15 MB RAM total!
// ============================================================================
#[derive(Clone, Copy, Debug)]
struct DD {
    hi: f64,
    lo: f64,
}

#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}

#[inline(always)]
fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let e = b - (s - a);
    (s, e)
}

#[inline(always)]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let e = a.mul_add(b, -p);
    (p, e)
}

impl DD {
    #[inline(always)]
    fn from_f64(x: f64) -> Self { Self { hi: x, lo: 0.0 } }
    #[inline(always)]
    fn from_i64(x: i64) -> Self {
        let hi = x as f64;
        let lo = (x - (hi as i64)) as f64;
        Self { hi, lo }
    }
    #[inline(always)]
    fn zero() -> Self { Self { hi: 0.0, lo: 0.0 } }
    #[inline(always)]
    fn one() -> Self { Self { hi: 1.0, lo: 0.0 } }
    #[inline(always)]
    fn neg(self) -> Self { Self { hi: -self.hi, lo: -self.lo } }
    #[inline(always)]
    fn abs(self) -> Self { if self.hi < 0.0 { self.neg() } else { self } }

    #[inline(always)]
    fn add(self, o: Self) -> Self {
        let (s1, s2) = two_sum(self.hi, o.hi);
        let (t1, t2) = two_sum(self.lo, o.lo);
        let (s1, s2) = quick_two_sum(s1, s2 + t1);
        let (hi, lo) = quick_two_sum(s1, s2 + t2);
        Self { hi, lo }
    }
    #[inline(always)]
    fn sub(self, o: Self) -> Self { self.add(o.neg()) }
    #[inline(always)]
    fn mul(self, o: Self) -> Self {
        let (p1, p2) = two_prod(self.hi, o.hi);
        let p2 = p2 + self.hi * o.lo + self.lo * o.hi;
        let (hi, lo) = quick_two_sum(p1, p2);
        Self { hi, lo }
    }
    #[inline(always)]
    fn div(self, o: Self) -> Self {
        let q1 = self.hi / o.hi;
        let r = self.sub(o.mul(DD::from_f64(q1)));
        let q2 = r.hi / o.hi;
        let r2 = r.sub(o.mul(DD::from_f64(q2)));
        let q3 = r2.hi / o.hi;
        let (s, e) = quick_two_sum(q1, q2);
        DD { hi: s, lo: e }.add(DD::from_f64(q3))
    }
    fn sqrt(self) -> Self {
        if self.hi <= 0.0 { return DD::zero(); }
        let x = DD::from_f64(self.hi.sqrt());
        // 1 Newton step doubles 53 bits -> 106 bits!
        x.add(self.div(x)).mul(DD::from_f64(0.5))
    }
    fn to_i128_scaled(self, scale: DD) -> Option<i128> {
        let v = self.mul(scale);
        if !v.hi.is_finite() || v.hi.abs() > 1e35 { return None; }
        // Split hi + lo cleanly into i128!
        let hi_i = v.hi.trunc() as i128;
        let rem = (v.hi - (hi_i as f64)) + v.lo;
        Some(hi_i + rem.round() as i128)
    }
}

fn parse_dd(s: &str) -> DD {
    // Parse decimal string with <= 34 digits into exact DD
    let neg = s.starts_with('-');
    let body = if neg { &s[1..] } else { s };
    let parts: Vec<&str> = body.split('.').collect();
    let int_part = parts[0];
    let frac_part = if parts.len() > 1 { parts[1] } else { "" };
    let mut val = DD::zero();
    let ten = DD::from_f64(10.0);
    for ch in int_part.chars() {
        let d = (ch as u8 - b'0') as f64;
        val = val.mul(ten).add(DD::from_f64(d));
    }
    let mut p = DD::one();
    for ch in frac_part.chars().take(34) {
        p = p.div(ten);
        let d = (ch as u8 - b'0') as f64;
        val = val.add(p.mul(DD::from_f64(d)));
    }
    if neg { val.neg() } else { val }
}

// ============================================================================
// Fast 5D LLL over i128 / f64 (detects short integer relations to < 1e-22!)
// ============================================================================
fn lll_find_rel(vals: &[DD], scale_pow10: u32, max_norm: f64, eps_hi: f64) -> Option<(Vec<i64>, f64)> {
    let n = vals.len();
    // Scale = 10^scale_pow10
    let mut scale = DD::one();
    for _ in 0..scale_pow10 { scale = scale.mul(DD::from_f64(10.0)); }
    let mut ints = Vec::with_capacity(n);
    for &v in vals {
        ints.push(v.to_i128_scaled(scale)?);
    }

    // Build basis B: n rows of length n+1 in f64 (and exact i128 copy for integer rows!)
    // Use i128 for basis vectors, f64 for Gram-Schmidt mu!
    let mut b: Vec<Vec<i128>> = vec![vec![0; n + 1]; n];
    for i in 0..n {
        b[i][i] = 1;
        b[i][n] = ints[i];
    }

    let mut bstar = vec![vec![0.0_f64; n + 1]; n];
    let mut mu = vec![vec![0.0_f64; n]; n];
    let mut bnorm = vec![0.0_f64; n];

    let gs = |b: &Vec<Vec<i128>>, bstar: &mut Vec<Vec<f64>>, mu: &mut Vec<Vec<f64>>, bnorm: &mut Vec<f64>| {
        for i in 0..n {
            for c in 0..=n { bstar[i][c] = b[i][c] as f64; }
            for j in 0..i {
                if bnorm[j] > 1e-30 {
                    let mut dot = 0.0;
                    for c in 0..=n { dot += (b[i][c] as f64) * bstar[j][c]; }
                    mu[i][j] = dot / bnorm[j];
                    for c in 0..=n { bstar[i][c] -= mu[i][j] * bstar[j][c]; }
                } else {
                    mu[i][j] = 0.0;
                }
            }
            let mut nm = 0.0;
            for c in 0..=n { nm += bstar[i][c] * bstar[i][c]; }
            bnorm[i] = nm;
        }
    };

    gs(&b, &mut bstar, &mut mu, &mut bnorm);
    let mut k = 1_usize;
    let mut iter = 0_usize;
    while k < n && iter < 400 {
        iter += 1;
        let mut changed = false;
        for j in (0..k).rev() {
            if mu[k][j].abs() > 0.500001 {
                let q = mu[k][j].round() as i128;
                if q != 0 {
                    for c in 0..=n {
                        b[k][c] = b[k][c].saturating_sub(q.saturating_mul(b[j][c]));
                    }
                    changed = true;
                }
            }
        }
        if changed { gs(&b, &mut bstar, &mut mu, &mut bnorm); }
        if bnorm[k] >= (0.75 - mu[k][k - 1] * mu[k][k - 1]) * bnorm[k - 1] {
            k += 1;
        } else {
            b.swap(k, k - 1);
            gs(&b, &mut bstar, &mut mu, &mut bnorm);
            k = if k > 1 { k - 1 } else { 1 };
        }
    }

    // Check every reduced row for a genuine DD residual < eps_hi with c[0] != 0!
    let mut best: Option<(Vec<i64>, f64)> = None;
    for row in &b {
        if row[0] == 0 { continue; }
        let mut norm2 = 0.0_f64;
        let mut ok = true;
        let mut coeffs = Vec::with_capacity(n);
        for i in 0..n {
            if row[i].abs() > 100_000 { ok = false; break; }
            let ci = row[i] as i64;
            coeffs.push(ci);
            norm2 += (ci as f64) * (ci as f64);
        }
        if !ok || norm2.sqrt() > max_norm { continue; }
        // Require at least one series term nonzero (indices >= 1 or >= 2)
        let mut sum = DD::zero();
        for i in 0..n {
            sum = sum.add(vals[i].mul(DD::from_i64(coeffs[i])));
        }
        let err = sum.hi.abs();
        if err < eps_hi {
            if coeffs[0] < 0 { for c in &mut coeffs { *c = -*c; } }
            if best.as_ref().map_or(true, |b| norm2 < b.0.iter().map(|&x| (x as f64).powi(2)).sum::<f64>()) {
                best = Some((coeffs, err));
            }
        }
    }
    best
}

// ============================================================================
// HYPERGEOMETRIC KERNEL SPECIFICATION
// Term: T_k = z^k / k^m * B1(k)^{e1} * B2(k)^{e2} * B3(k)^{e3} * B4(k)^{e4}
//   B1(k) = (1/2)_k / (1)_k           = C(2k,k) / 4^k
//   B2(k) = (1/3)_k(2/3)_k / (1)_k^2  = C(3k,k)C(2k,k) / 27^k
//   B3(k) = (1/4)_k(3/4)_k / (1)_k^2  = C(4k,2k)C(2k,k) / 64^k
//   B4(k) = (1/6)_k(5/6)_k / (1)_k^2  = C(6k,3k)C(3k,k) / 432^k
//   B5(k) = (1/2)_k^2 / (1)_k^2       (already covered by e1 += 2!)
// ============================================================================

#[derive(Clone, Debug)]
struct BaseDef {
    name: &'static str,
    val: DD,
}

#[derive(Clone, Debug)]
struct Hit {
    target: &'static str,
    e: [i8; 4],
    m: i8,
    base_name: &'static str,
    coeffs: Vec<i64>,
    err: f64,
    terms_used: usize,
}

fn eval_moments_dd(e: [i8; 4], m: i8, z: DD) -> Option<([DD; 4], usize)> {
    // Evaluate S_r = sum_{k=1}^K k^r * (z^k / k^m) * prod B_i(k)^{e_i} for r=0,1,2,3
    // Converges fast if |z| < 0.85! Stop when |term * k^3| < 1e-30 (max 320 terms).
    let mut s = [DD::zero(); 4];
    let mut prod = DD::one(); // holds z^k * prod B_i(k)^{e_i}
    let eps = 1e-29;

    // Precompute recurrence step factor for each B_i:
    // B1(k)/B1(k-1) = (2k-1)/(2k)
    // B2(k)/B2(k-1) = (3k-2)(3k-1)/(9k^2)
    // B3(k)/B3(k-1) = (4k-3)(4k-1)/(16k^2)
    // B4(k)/B4(k-1) = (6k-5)(6k-1)/(36k^2)
    let mut converged = false;
    let mut used = 0_usize;
    for k in 1..=320_usize {
        let kf = k as f64;
        let k_dd = DD::from_f64(kf);
        let k2_dd = DD::from_f64(kf * kf);
        prod = prod.mul(z);

        // e[0]: B1
        if e[0] != 0 {
            let num = DD::from_f64(2.0 * kf - 1.0);
            let den = DD::from_f64(2.0 * kf);
            let r = num.div(den);
            let cnt = e[0].unsigned_abs();
            for _ in 0..cnt { prod = if e[0] > 0 { prod.mul(r) } else { prod.div(r) }; }
        }
        // e[1]: B2
        if e[1] != 0 {
            let num = DD::from_f64((3.0 * kf - 2.0) * (3.0 * kf - 1.0));
            let den = k2_dd.mul(DD::from_f64(9.0));
            let r = num.div(den);
            let cnt = e[1].unsigned_abs();
            for _ in 0..cnt { prod = if e[1] > 0 { prod.mul(r) } else { prod.div(r) }; }
        }
        // e[2]: B3
        if e[2] != 0 {
            let num = DD::from_f64((4.0 * kf - 3.0) * (4.0 * kf - 1.0));
            let den = k2_dd.mul(DD::from_f64(16.0));
            let r = num.div(den);
            let cnt = e[2].unsigned_abs();
            for _ in 0..cnt { prod = if e[2] > 0 { prod.mul(r) } else { prod.div(r) }; }
        }
        // e[3]: B4
        if e[3] != 0 {
            let num = DD::from_f64((6.0 * kf - 5.0) * (6.0 * kf - 1.0));
            let den = k2_dd.mul(DD::from_f64(36.0));
            let r = num.div(den);
            let cnt = e[3].unsigned_abs();
            for _ in 0..cnt { prod = if e[3] > 0 { prod.mul(r) } else { prod.div(r) }; }
        }

        if !prod.hi.is_finite() || prod.hi.abs() > 1e100 { return None; }
        // Divide by k^m
        let mut km = DD::one();
        for _ in 0..m { km = km.mul(k_dd); }
        let t0 = prod.div(km);
        let t1 = t0.mul(k_dd);
        let t2 = t1.mul(k_dd);
        let t3 = t2.mul(k_dd);
        s[0] = s[0].add(t0);
        s[1] = s[1].add(t1);
        s[2] = s[2].add(t2);
        s[3] = s[3].add(t3);
        used = k;
        if k >= 12 && t3.hi.abs() < eps * (s[3].hi.abs().max(1.0)) {
            converged = true;
            break;
        }
    }
    if !converged { return None; }
    // Avoid trivial near-zero sums
    if s[0].hi.abs() < 1e-18 { return None; }
    Some((s, used))
}

fn main() {
    let t0 = Instant::now();
    // 34-digit constants
    let z3 = parse_dd("1.2020569031595942853997381615114500");
    let pi = parse_dd("3.1415926535897932384626433832795029");
    let pi3 = pi.mul(pi).mul(pi);
    let z3_over_pi3 = z3.div(pi3);
    let inv_z3 = DD::one().div(z3);
    let phi = DD::one().add(DD::from_f64(5.0).sqrt()).mul(DD::from_f64(0.5));
    let inv_phi = DD::one().div(phi);
    let silver = DD::one().add(DD::from_f64(2.0).sqrt());
    let inv_sil = DD::one().div(silver);
    let rho = parse_dd("1.3247179572447460259609088544780973");
    let inv_rho = DD::one().div(rho);

    // Rich catalog of rational + algebraic bases |z| <= 0.78 (fast 31-digit convergence in <= 320 terms!)
    let mut bases: Vec<BaseDef> = Vec::new();
    let rat_bases: &[(&str, i64, i64)] = &[
        ("1/2", 1, 2), ("-1/2", -1, 2),
        ("1/3", 1, 3), ("-1/3", -1, 3),
        ("1/4", 1, 4), ("-1/4", -1, 4),
        ("1/5", 1, 5), ("-1/5", -1, 5),
        ("1/6", 1, 6), ("-1/6", -1, 6),
        ("1/8", 1, 8), ("-1/8", -1, 8),
        ("1/9", 1, 9), ("-1/9", -1, 9),
        ("1/12", 1, 12), ("-1/12", -1, 12),
        ("1/16", 1, 16), ("-1/16", -1, 16),
        ("1/24", 1, 24), ("-1/24", -1, 24),
        ("1/27", 1, 27), ("-1/27", -1, 27),
        ("1/32", 1, 32), ("-1/32", -1, 32),
        ("1/48", 1, 48), ("-1/48", -1, 48),
        ("1/64", 1, 64), ("-1/64", -1, 64),
        ("1/81", 1, 81), ("-1/81", -1, 81),
        ("1/108", 1, 108), ("-1/108", -1, 108),
        ("1/128", 1, 128), ("-1/128", -1, 128),
        ("1/256", 1, 256), ("-1/256", -1, 256),
        ("1/432", 1, 432), ("-1/432", -1, 432),
        ("1/1024", 1, 1024), ("-1/1024", -1, 1024),
        ("1/1728", 1, 1728), ("-1/1728", -1, 1728),
        ("2/3", 2, 3), ("-2/3", -2, 3),
        ("3/4", 3, 4), ("-3/4", -3, 4),
        ("2/9", 2, 9), ("-2/9", -2, 9),
        ("4/9", 4, 9), ("-4/9", -4, 9),
        ("4/27", 4, 27), ("-4/27", -4, 27),
        ("8/27", 8, 27), ("-8/27", -8, 27),
        ("16/27", 16, 27), ("-16/27", -16, 27),
        ("27/64", 27, 64), ("-27/64", -27, 64),
        ("9/16", 9, 16), ("-9/16", -9, 16),
        ("1/7", 1, 7), ("-1/7", -1, 7),
        ("1/10", 1, 10), ("-1/10", -1, 10),
        ("1/18", 1, 18), ("-1/18", -1, 18),
        ("1/36", 1, 36), ("-1/36", -1, 36),
        ("1/72", 1, 72), ("-1/72", -1, 72),
        ("1/144", 1, 144), ("-1/144", -1, 144),
        ("1/512", 1, 512), ("-1/512", -1, 512),
        ("1/2048", 1, 2048), ("-1/2048", -1, 2048),
        ("1/4096", 1, 4096), ("-1/4096", -1, 4096),
    ];
    for &(nm, p, q) in rat_bases {
        bases.push(BaseDef { name: nm, val: DD::from_i64(p).div(DD::from_i64(q)) });
    }
    // Algebraic bases: golden, silver, plastic!
    let p1 = inv_phi;
    let p2 = p1.mul(p1);
    let p3 = p2.mul(p1);
    let p4 = p2.mul(p2);
    let p5 = p4.mul(p1);
    let p6 = p3.mul(p3);
    let s2 = inv_sil.mul(inv_sil);
    let s4 = s2.mul(s2);
    let r3 = inv_rho.mul(inv_rho).mul(inv_rho);
    let r4 = r3.mul(inv_rho);
    let r5 = r4.mul(inv_rho);
    for (nm, val) in [
        ("phi^-1", p1), ("-phi^-1", p1.neg()),
        ("phi^-2", p2), ("-phi^-2", p2.neg()),
        ("phi^-3", p3), ("-phi^-3", p3.neg()),
        ("phi^-4", p4), ("-phi^-4", p4.neg()),
        ("phi^-5", p5), ("-phi^-5", p5.neg()),
        ("phi^-6", p6), ("-phi^-6", p6.neg()),
        ("sil^-2", s2), ("-sil^-2", s2.neg()),
        ("sil^-4", s4), ("-sil^-4", s4.neg()),
        ("rho^-3", r3), ("-rho^-3", r3.neg()),
        ("rho^-4", r4), ("-rho^-4", r4.neg()),
        ("rho^-5", r5), ("-rho^-5", r5.neg()),
    ] {
        bases.push(BaseDef { name: nm, val });
    }

    // Generate exponent tuples e = [e1, e2, e3, e4] covering rich fractional Pochhammer combos!
    let mut exps: Vec<[i8; 4]> = Vec::new();
    for e1 in -5_i8..=5 {
        for e2 in -3_i8..=3 {
            for e3 in -3_i8..=3 {
                for e4 in -2_i8..=2 {
                    let l1 = e1.abs() + e2.abs() + e3.abs() + e4.abs();
                    if l1 > 6 { continue; }
                    // Avoid overly mixed 4-factor combos unless small
                    let nz = (e1 != 0) as i8 + (e2 != 0) as i8 + (e3 != 0) as i8 + (e4 != 0) as i8;
                    if nz > 2 && l1 > 4 { continue; }
                    exps.push([e1, e2, e3, e4]);
                }
            }
        }
    }

    let mut jobs: Vec<([i8; 4], i8, usize)> = Vec::new();
    for &e in &exps {
        for m in 3_i8..=5 {
            for bi in 0..bases.len() {
                jobs.push((e, m, bi));
            }
        }
    }
    println!("=== 16-CORE RUST FRACTIONAL-POCHHAMMER HYPERGEOMETRIC HUNTER ===");
    println!("Exponent tuples: {} | Bases: {} | Total (e, m, z) kernels: {}",
             exps.len(), bases.len(), jobs.len());

    let jobs_arc = Arc::new(jobs);
    let bases_arc = Arc::new(bases);
    let hits_m = Arc::new(Mutex::new(Vec::<Hit>::new()));
    let n_threads = 16;
    let chunk = (jobs_arc.len() + n_threads - 1) / n_threads;

    let mut handles = Vec::new();
    for tid in 0..n_threads {
        let j_ref = Arc::clone(&jobs_arc);
        let b_ref = Arc::clone(&bases_arc);
        let h_ref = Arc::clone(&hits_m);
        handles.push(thread::spawn(move || {
            let lo = tid * chunk;
            let hi = (lo + chunk).min(j_ref.len());
            let mut local_hits = Vec::new();
            for idx in lo..hi {
                let (e, m, bi) = j_ref[idx];
                let base = &b_ref[bi];
                let Some((s, used)) = eval_moments_dd(e, m, base.val) else { continue; };

                // Target 1: PURE zeta(3) = sum (c1 + c2*k + c3*k^2 + c4*k^3) * T_k
                if let Some((coeffs, err)) = lll_find_rel(&[z3, s[0], s[1], s[2], s[3]], 25, 2500.0, 1e-21) {
                    // Reject trivial if all series coeffs are 0
                    if coeffs[1..].iter().any(|&c| c != 0) {
                        local_hits.push(Hit { target: "zeta(3)", e, m, base_name: base.name, coeffs, err, terms_used: used });
                    }
                }
                // Target 2: PURE pi^3 = sum (c1 + c2*k + c3*k^2 + c4*k^3) * T_k
                if let Some((coeffs, err)) = lll_find_rel(&[pi3, s[0], s[1], s[2], s[3]], 25, 2500.0, 1e-21) {
                    if coeffs[1..].iter().any(|&c| c != 0) {
                        local_hits.push(Hit { target: "pi^3", e, m, base_name: base.name, coeffs, err, terms_used: used });
                    }
                }
                // Target 3: zeta(3)/pi^3 = sum (c1 + c2*k + c3*k^2) * T_k  (plus rational const!)
                if let Some((coeffs, err)) = lll_find_rel(&[z3_over_pi3, DD::one(), s[0], s[1], s[2]], 25, 600.0, 1e-21) {
                    if coeffs[2..].iter().any(|&c| c != 0) {
                        local_hits.push(Hit { target: "zeta(3)/pi^3", e, m, base_name: base.name, coeffs, err, terms_used: used });
                    }
                }
                // Target 4: 1/zeta(3) = r + sum (c1 + c2*k + c3*k^2) * T_k
                if let Some((coeffs, err)) = lll_find_rel(&[inv_z3, DD::one(), s[0], s[1], s[2]], 25, 600.0, 1e-21) {
                    if coeffs[2..].iter().any(|&c| c != 0) {
                        local_hits.push(Hit { target: "1/zeta(3)", e, m, base_name: base.name, coeffs, err, terms_used: used });
                    }
                }
            }
            h_ref.lock().unwrap().extend(local_hits);
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut hits = Arc::try_unwrap(hits_m).unwrap().into_inner().unwrap();
    hits.sort_by(|a, b| (a.target, a.e, a.m, a.base_name).cmp(&(b.target, b.e, b.m, b.base_name)));
    println!("\nScanned {} kernels in {:.2}s! Found {} candidate hits (err < 1e-21):",
             jobs_arc.len(), t0.elapsed().as_secs_f64(), hits.len());
    for h in &hits {
        println!("  HIT [{:12}] e={:?} m={} z={:>8} terms={:3} err={:.2e} coeffs={:?}",
                 h.target, h.e, h.m, h.base_name, h.terms_used, h.err, h.coeffs);
    }
}
