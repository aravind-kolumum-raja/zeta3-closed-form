use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

#[derive(Clone, Copy, Debug)]
struct DD { hi: f64, lo: f64 }

#[inline(always)] fn two_sum(a: f64, b: f64) -> (f64, f64) { let s = a + b; let v = s - a; (s, (a - (s - v)) + (b - v)) }
#[inline(always)] fn quick_two_sum(a: f64, b: f64) -> (f64, f64) { let s = a + b; (s, b - (s - a)) }
#[inline(always)] fn two_prod(a: f64, b: f64) -> (f64, f64) { let p = a * b; (p, a.mul_add(b, -p)) }

impl DD {
    #[inline(always)] fn from_f64(x: f64) -> Self { Self { hi: x, lo: 0.0 } }
    #[inline(always)] fn from_i64(x: i64) -> Self { let hi = x as f64; Self { hi, lo: (x - (hi as i64)) as f64 } }
    #[inline(always)] fn zero() -> Self { Self { hi: 0.0, lo: 0.0 } }
    #[inline(always)] fn one() -> Self { Self { hi: 1.0, lo: 0.0 } }
    #[inline(always)] fn neg(self) -> Self { Self { hi: -self.hi, lo: -self.lo } }
    #[inline(always)] fn add(self, o: Self) -> Self {
        let (s1, s2) = two_sum(self.hi, o.hi);
        let (t1, t2) = two_sum(self.lo, o.lo);
        let (s1, s2) = quick_two_sum(s1, s2 + t1);
        let (hi, lo) = quick_two_sum(s1, s2 + t2);
        Self { hi, lo }
    }
    #[inline(always)] fn sub(self, o: Self) -> Self { self.add(o.neg()) }
    #[inline(always)] fn mul(self, o: Self) -> Self {
        let (p1, p2) = two_prod(self.hi, o.hi);
        let (hi, lo) = quick_two_sum(p1, p2 + self.hi * o.lo + self.lo * o.hi);
        Self { hi, lo }
    }
    #[inline(always)] fn div(self, o: Self) -> Self {
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
        x.add(self.div(x)).mul(DD::from_f64(0.5))
    }
    fn to_i128_scaled(self, scale: DD) -> Option<i128> {
        let v = self.mul(scale);
        if !v.hi.is_finite() || v.hi.abs() > 1e34 { return None; }
        let hi_i = v.hi.trunc() as i128;
        let rem = (v.hi - (hi_i as f64)) + v.lo;
        Some(hi_i + rem.round() as i128)
    }
}

fn parse_dd(s: &str) -> DD {
    let neg = s.starts_with('-');
    let body = if neg { &s[1..] } else { s };
    let parts: Vec<&str> = body.split('.').collect();
    let mut val = DD::zero();
    let ten = DD::from_f64(10.0);
    for ch in parts[0].chars() { val = val.mul(ten).add(DD::from_f64((ch as u8 - b'0') as f64)); }
    if parts.len() > 1 {
        let mut p = DD::one();
        for ch in parts[1].chars().take(34) {
            p = p.div(ten);
            val = val.add(p.mul(DD::from_f64((ch as u8 - b'0') as f64)));
        }
    }
    if neg { val.neg() } else { val }
}

fn lll_find_rel(vals: &[DD], scale_pow10: u32, max_norm: f64, eps_hi: f64, series_start: usize) -> Option<(Vec<i64>, f64)> {
    let n = vals.len();
    let mut scale = DD::one();
    for _ in 0..scale_pow10 { scale = scale.mul(DD::from_f64(10.0)); }
    let mut ints = Vec::with_capacity(n);
    for &v in vals { ints.push(v.to_i128_scaled(scale)?); }

    let mut b: Vec<Vec<i128>> = vec![vec![0; n + 1]; n];
    for i in 0..n { b[i][i] = 1; b[i][n] = ints[i]; }
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
                } else { mu[i][j] = 0.0; }
            }
            let mut nm = 0.0;
            for c in 0..=n { nm += bstar[i][c] * bstar[i][c]; }
            bnorm[i] = nm;
        }
    };

    gs(&b, &mut bstar, &mut mu, &mut bnorm);
    let mut k = 1_usize; let mut iter = 0_usize;
    while k < n && iter < 350 {
        iter += 1;
        let mut changed = false;
        for j in (0..k).rev() {
            if mu[k][j].abs() > 0.500001 {
                let q = mu[k][j].round() as i128;
                if q != 0 {
                    for c in 0..=n { b[k][c] = b[k][c].saturating_sub(q.saturating_mul(b[j][c])); }
                    changed = true;
                }
            }
        }
        if changed { gs(&b, &mut bstar, &mut mu, &mut bnorm); }
        if bnorm[k] >= (0.75 - mu[k][k - 1] * mu[k][k - 1]) * bnorm[k - 1] { k += 1; }
        else { b.swap(k, k - 1); gs(&b, &mut bstar, &mut mu, &mut bnorm); k = if k > 1 { k - 1 } else { 1 }; }
    }

    let mut best: Option<(Vec<i64>, f64)> = None;
    for row in &b {
        if row[0] == 0 { continue; }
        let mut norm2 = 0.0_f64; let mut ok = true;
        let mut coeffs = Vec::with_capacity(n);
        for i in 0..n {
            if row[i].abs() > 100_000 { ok = false; break; }
            let ci = row[i] as i64;
            coeffs.push(ci);
            norm2 += (ci as f64) * (ci as f64);
        }
        if !ok || norm2.sqrt() > max_norm { continue; }
        if !coeffs[series_start..].iter().any(|&c| c != 0) { continue; }
        let mut sum = DD::zero();
        for i in 0..n { sum = sum.add(vals[i].mul(DD::from_i64(coeffs[i]))); }
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
// Denominator & Pochhammer families:
//   DenomKind:
//     0: k^3       (k >= 1)
//     1: k^4       (k >= 1)
//     2: k^5       (k >= 1)
//     3: (2k+1)^3  (k >= 0) -- odd cubic! (gives pi^3 & Sun/Guillera zeta(3)!)
//     4: (2k+1)^5  (k >= 0) -- odd quintic!
//     5: (2k-1)^3  (k >= 1) -- shifted odd cubic!
//     6: k^3*(2k-1)(2k+1)? No, k^2*(2k-1)^2 (k >= 1) and k*(2k-1)^2 (k >= 1)
//     7: 1         (k >= 0) -- pure Ramanujan-Sato (Guillera 1/pi^2, 1/zeta(3), etc.!)
//   Pochhammer ratios (all weight-0 normalized so convergence radius is |z| < 1!):
//     F0(k) = (1/2)_k / (1)_k            = C(2k,k) / 4^k
//     F1(k) = (1/3)_k(2/3)_k / (1)_k^2   = C(3k,k)C(2k,k) / 27^k
//     F2(k) = (1/4)_k(3/4)_k / (1)_k^2   = C(4k,2k)C(2k,k) / 64^k
//     F3(k) = (1/6)_k(5/6)_k / (1)_k^2   = C(6k,3k)C(3k,k) / 432^k
//     Note: C(3k,k)/(27/4)^k = F1/F0 (already covered by e=[e0, e1, ...])!
//     Note: C(4k,2k)/16^k    = F2/F0 (already covered by e=[e0, ..., e2])!
// ============================================================================

fn denom_name(dk: u8) -> &'static str {
    match dk {
        0 => "k^3",
        1 => "k^4",
        2 => "k^5",
        3 => "(2k+1)^3[k>=0]",
        4 => "(2k+1)^5[k>=0]",
        5 => "(2k-1)^3[k>=1]",
        6 => "k^2(2k-1)[k>=1]",
        7 => "1[k>=0]",
        _ => "?",
    }
}

fn eval_moments(e: [i8; 4], dk: u8, z: DD) -> Option<([DD; 4], usize)> {
    let k_start = match dk { 3 | 4 | 7 => 0_usize, _ => 1_usize };
    let mut s = [DD::zero(); 4];
    let mut prod = DD::one();
    let eps = 1e-29;
    let mut converged = false;
    let mut used = 0_usize;

    for k in 0..=260_usize {
        if k > 0 {
            let kf = k as f64;
            let k2 = DD::from_f64(kf * kf);
            prod = prod.mul(z);
            if e[0] != 0 {
                let r = DD::from_f64(2.0 * kf - 1.0).div(DD::from_f64(2.0 * kf));
                for _ in 0..e[0].unsigned_abs() { prod = if e[0] > 0 { prod.mul(r) } else { prod.div(r) }; }
            }
            if e[1] != 0 {
                let r = DD::from_f64((3.0 * kf - 2.0) * (3.0 * kf - 1.0)).div(k2.mul(DD::from_f64(9.0)));
                for _ in 0..e[1].unsigned_abs() { prod = if e[1] > 0 { prod.mul(r) } else { prod.div(r) }; }
            }
            if e[2] != 0 {
                let r = DD::from_f64((4.0 * kf - 3.0) * (4.0 * kf - 1.0)).div(k2.mul(DD::from_f64(16.0)));
                for _ in 0..e[2].unsigned_abs() { prod = if e[2] > 0 { prod.mul(r) } else { prod.div(r) }; }
            }
            if e[3] != 0 {
                let r = DD::from_f64((6.0 * kf - 5.0) * (6.0 * kf - 1.0)).div(k2.mul(DD::from_f64(36.0)));
                for _ in 0..e[3].unsigned_abs() { prod = if e[3] > 0 { prod.mul(r) } else { prod.div(r) }; }
            }
            if !prod.hi.is_finite() || prod.hi.abs() > 1e80 { return None; }
        }
        if k < k_start { continue; }

        let kf = k as f64;
        let k_dd = DD::from_f64(kf);
        let den = match dk {
            0 => { let d = k_dd.mul(k_dd).mul(k_dd); d }
            1 => { let d2 = k_dd.mul(k_dd); d2.mul(d2) }
            2 => { let d2 = k_dd.mul(k_dd); d2.mul(d2).mul(k_dd) }
            3 => { let o = DD::from_f64(2.0 * kf + 1.0); o.mul(o).mul(o) }
            4 => { let o = DD::from_f64(2.0 * kf + 1.0); let o2 = o.mul(o); o2.mul(o2).mul(o) }
            5 => { let o = DD::from_f64(2.0 * kf - 1.0); o.mul(o).mul(o) }
            6 => { k_dd.mul(k_dd).mul(DD::from_f64(2.0 * kf - 1.0)) }
            7 => DD::one(),
            _ => unreachable!(),
        };
        let t0 = prod.div(den);
        let t1 = t0.mul(k_dd);
        let t2 = t1.mul(k_dd);
        let t3 = t2.mul(k_dd);
        s[0] = s[0].add(t0);
        s[1] = s[1].add(t1);
        s[2] = s[2].add(t2);
        s[3] = s[3].add(t3);
        used = k;
        if k >= 12 && t3.hi.abs() < eps * (s[3].hi.abs().max(1.0)) && t0.hi.abs() < eps * (s[0].hi.abs().max(1.0)) {
            converged = true;
            break;
        }
    }
    if !converged || s[0].hi.abs() < 1e-18 { return None; }
    Some((s, used))
}

#[derive(Clone, Debug)]
struct Hit {
    target: &'static str,
    e: [i8; 4],
    dk: u8,
    base: &'static str,
    coeffs: Vec<i64>,
    err: f64,
    terms: usize,
}

fn main() {
    let t0 = Instant::now();
    let z3 = parse_dd("1.2020569031595942853997381615114500");
    let pi = parse_dd("3.1415926535897932384626433832795029");
    let pi2 = pi.mul(pi);
    let pi3 = pi2.mul(pi);
    let inv_pi2 = DD::one().div(pi2);
    let ln2 = parse_dd("0.6931471805599453094172321214581766");
    let pi2_ln2 = pi2.mul(ln2);
    let z3_over_pi3 = z3.div(pi3);
    let phi = DD::one().add(DD::from_f64(5.0).sqrt()).mul(DD::from_f64(0.5));
    let inv_phi = DD::one().div(phi);
    let rho = parse_dd("1.3247179572447460259609088544780973");
    let inv_rho = DD::one().div(rho);

    let mut bases: Vec<(&'static str, DD)> = Vec::new();
    let rats: &[(&str, i64, i64)] = &[
        ("1/2", 1, 2), ("-1/2", -1, 2), ("1/3", 1, 3), ("-1/3", -1, 3),
        ("1/4", 1, 4), ("-1/4", -1, 4), ("1/8", 1, 8), ("-1/8", -1, 8),
        ("1/9", 1, 9), ("-1/9", -1, 9), ("1/16", 1, 16), ("-1/16", -1, 16),
        ("1/27", 1, 27), ("-1/27", -1, 27), ("1/32", 1, 32), ("-1/32", -1, 32),
        ("1/64", 1, 64), ("-1/64", -1, 64), ("1/81", 1, 81), ("-1/81", -1, 81),
        ("1/128", 1, 128), ("-1/128", -1, 128), ("1/256", 1, 256), ("-1/256", -1, 256),
        ("1/512", 1, 512), ("-1/512", -1, 512), ("1/1024", 1, 1024), ("-1/1024", -1, 1024),
        ("1/4096", 1, 4096), ("-1/4096", -1, 4096),
        ("1/12", 1, 12), ("-1/12", -1, 12), ("1/24", 1, 24), ("-1/24", -1, 24),
        ("1/48", 1, 48), ("-1/48", -1, 48), ("1/108", 1, 108), ("-1/108", -1, 108),
        ("1/432", 1, 432), ("-1/432", -1, 432), ("1/1728", 1, 1728), ("-1/1728", -1, 1728),
        ("4/27", 4, 27), ("-4/27", -4, 27), ("27/64", 27, 64), ("-27/64", -27, 64),
        ("16/27", 16, 27), ("-16/27", -16, 27), ("9/16", 9, 16), ("-9/16", -9, 16),
        ("1/5", 1, 5), ("-1/5", -1, 5), ("1/25", 1, 25), ("-1/25", -1, 25),
        ("1/6", 1, 6), ("-1/6", -1, 6), ("1/36", 1, 36), ("-1/36", -1, 36),
        ("1/7", 1, 7), ("-1/7", -1, 7), ("1/49", 1, 49), ("-1/49", -1, 49),
        ("3/4", 3, 4), ("-3/4", -3, 4), ("2/3", 2, 3), ("-2/3", -2, 3),
        ("1/18", 1, 18), ("-1/18", -1, 18), ("1/72", 1, 72), ("-1/72", -1, 72),
    ];
    for &(nm, p, q) in rats { bases.push((nm, DD::from_i64(p).div(DD::from_i64(q)))); }
    let p1 = inv_phi; let p2 = p1.mul(p1); let p3 = p2.mul(p1); let p4 = p2.mul(p2);
    let r3 = inv_rho.mul(inv_rho).mul(inv_rho); let r4 = r3.mul(inv_rho);
    for (nm, v) in [
        ("phi^-1", p1), ("-phi^-1", p1.neg()), ("phi^-2", p2), ("-phi^-2", p2.neg()),
        ("phi^-3", p3), ("-phi^-3", p3.neg()), ("phi^-4", p4), ("-phi^-4", p4.neg()),
        ("rho^-3", r3), ("-rho^-3", r3.neg()), ("rho^-4", r4), ("-rho^-4", r4.neg()),
    ] { bases.push((nm, v)); }

    let mut exps: Vec<[i8; 4]> = Vec::new();
    for e0 in -6_i8..=6 {
        for e1 in -3_i8..=3 {
            for e2 in -3_i8..=3 {
                for e3 in -2_i8..=2 {
                    let l1 = e0.abs() + e1.abs() + e2.abs() + e3.abs();
                    if l1 > 6 { continue; }
                    let nz = (e0 != 0) as i8 + (e1 != 0) as i8 + (e2 != 0) as i8 + (e3 != 0) as i8;
                    if nz > 2 && l1 > 4 { continue; }
                    exps.push([e0, e1, e2, e3]);
                }
            }
        }
    }

    let mut jobs: Vec<([i8; 4], u8, usize)> = Vec::new();
    for &e in &exps {
        for dk in 0_u8..=7 {
            for bi in 0..bases.len() { jobs.push((e, dk, bi)); }
        }
    }
    println!("=== WZ v2: {} kernels (8 denom families x {} Pochhammer combos x {} bases) ===",
             jobs.len(), exps.len(), bases.len());

    let jobs_arc = Arc::new(jobs);
    let bases_arc = Arc::new(bases);
    let hits_m = Arc::new(Mutex::new(Vec::<Hit>::new()));
    let n_threads = 16;
    let chunk = (jobs_arc.len() + n_threads - 1) / n_threads;
    let mut handles = Vec::new();

    for tid in 0..n_threads {
        let jr = Arc::clone(&jobs_arc); let br = Arc::clone(&bases_arc); let hr = Arc::clone(&hits_m);
        handles.push(thread::spawn(move || {
            let lo = tid * chunk; let hi = (lo + chunk).min(jr.len());
            let mut local = Vec::new();
            for idx in lo..hi {
                let (e, dk, bi) = jr[idx];
                let (bname, z) = br[bi];
                let Some((s, used)) = eval_moments(e, dk, z) else { continue; };

                // 1. Pure zeta(3)
                if let Some((c, err)) = lll_find_rel(&[z3, s[0], s[1], s[2], s[3]], 25, 2500.0, 1e-21, 1) {
                    local.push(Hit { target: "zeta(3)", e, dk, base: bname, coeffs: c, err, terms: used });
                }
                // 2. Pure pi^3
                if let Some((c, err)) = lll_find_rel(&[pi3, s[0], s[1], s[2], s[3]], 25, 2500.0, 1e-21, 1) {
                    local.push(Hit { target: "pi^3", e, dk, base: bname, coeffs: c, err, terms: used });
                }
                // 3. Mixed zeta(3) + r*pi^2*ln(2) (only check if not already pure!)
                if let Some((c, err)) = lll_find_rel(&[z3, pi2_ln2, s[0], s[1], s[2]], 25, 800.0, 1e-21, 2) {
                    if c[1] != 0 {
                        local.push(Hit { target: "z3+pi2ln2", e, dk, base: bname, coeffs: c, err, terms: used });
                    }
                }
                // 4. Guillera 1/pi^2 benchmark (for dk == 7 pure Ramanujan-Sato!)
                if dk == 7 {
                    if let Some((c, err)) = lll_find_rel(&[inv_pi2, s[0], s[1], s[2], s[3]], 25, 2500.0, 1e-21, 1) {
                        local.push(Hit { target: "1/pi^2", e, dk, base: bname, coeffs: c, err, terms: used });
                    }
                }
                // 5. zeta(3)/pi^3
                if let Some((c, err)) = lll_find_rel(&[z3_over_pi3, DD::one(), s[0], s[1], s[2]], 25, 500.0, 1e-21, 2) {
                    local.push(Hit { target: "z3/pi^3", e, dk, base: bname, coeffs: c, err, terms: used });
                }
            }
            hr.lock().unwrap().extend(local);
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut hits = Arc::try_unwrap(hits_m).unwrap().into_inner().unwrap();
    hits.sort_by(|a, b| (a.target, a.dk, a.e, a.base).cmp(&(b.target, b.dk, b.e, b.base)));
    println!("\nScanned {} kernels in {:.2}s! Found {} candidate hits (err < 1e-21):",
             jobs_arc.len(), t0.elapsed().as_secs_f64(), hits.len());
    for h in &hits {
        println!("  HIT [{:9}] dk={:<15} e={:?} z={:>7} terms={:3} err={:.2e} coeffs={:?}",
                 h.target, denom_name(h.dk), h.e, h.base, h.terms, h.err, h.coeffs);
    }
}
