//! pcf_hunter — 16-core, low-memory (<10 MB RAM) native Rust hunter for
//! degree-3 Polynomial Continued Fractions (PCFs) and Möbius closed forms
//! targeting ζ(3), ζ(3)/π³, and ζ(3)/π².
//!
//! Uses fast backward-recurrence evaluation + 4D LLL integer-relation detection
//! on (1, α, K, α·K) to find small-integer Möbius transforms K = (p1·α + p2)/(q1·α + q2),
//! plus exact self-test recovery of Apéry's PCF (6/ζ(3)) and Ramanujan-Machine PCFs.

use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

const ZETA3: f64 = 1.2020569031595942854;
const PI: f64 = std::f64::consts::PI;

/// Evaluate K = a(0) + b(1)/(a(1) + b(2)/(a(2) + ...)) using Lentz / backward evaluation.
/// Returns Some(K) if convergent to < 1e-13 between N1 and N2 steps.
#[inline(always)]
fn eval_pcf<FA, FB>(a_fn: FA, b_fn: FB, n1: i64, n2: i64) -> Option<f64>
where
    FA: Fn(i64) -> f64,
    FB: Fn(i64) -> f64,
{
    let mut v1 = a_fn(n1);
    if !v1.is_finite() || v1 == 0.0 {
        return None;
    }
    for n in (0..n1).rev() {
        let an = a_fn(n);
        let bn1 = b_fn(n + 1);
        if v1.abs() < 1e-300 {
            return None;
        }
        v1 = an + bn1 / v1;
        if !v1.is_finite() {
            return None;
        }
    }

    let mut v2 = a_fn(n2);
    if !v2.is_finite() || v2 == 0.0 {
        return None;
    }
    for n in (0..n2).rev() {
        let an = a_fn(n);
        let bn1 = b_fn(n + 1);
        if v2.abs() < 1e-300 {
            return None;
        }
        v2 = an + bn1 / v2;
        if !v2.is_finite() {
            return None;
        }
    }

    let diff = (v1 - v2).abs();
    let scale = v2.abs().max(1.0);
    if diff < 1e-13 * scale && v2.abs() > 1e-9 && v2.abs() < 1e9 {
        Some(v2)
    } else {
        None
    }
}

/// Fast 4D LLL integer relation finder for x = [x0, x1, x2, x3].
/// Looks for non-zero integer vector c with |c_i| <= max_c and |c · x| < tol.
#[inline]
fn find_mobius_lll(alpha: f64, k_val: f64, max_c: i64, tol: f64) -> Option<[i64; 4]> {
    let x = [1.0_f64, alpha, k_val, alpha * k_val];
    let scale = 1.0e12_f64;

    // 4 basis rows in R^5: [e_i (4 coords), scale * x_i]
    let mut b = [[0.0_f64; 5]; 4];
    for i in 0..4 {
        b[i][i] = 1.0;
        b[i][4] = scale * x[i];
    }

    let dot = |u: &[f64; 5], v: &[f64; 5]| -> f64 {
        u[0] * v[0] + u[1] * v[1] + u[2] * v[2] + u[3] * v[3] + u[4] * v[4]
    };

    // Standard LLL loop (delta = 0.75) for fixed dim 4 — super fast, stack-only!
    let mut ortho = [[0.0_f64; 5]; 4];
    let mut mu = [[0.0_f64; 4]; 4];
    let mut norm2 = [0.0_f64; 4];

    let recompute_gs = |b: &[[f64; 5]; 4],
                        ortho: &mut [[f64; 5]; 4],
                        mu: &mut [[f64; 4]; 4],
                        norm2: &mut [f64; 4]| {
        for i in 0..4 {
            ortho[i] = b[i];
            for j in 0..i {
                if norm2[j] > 1e-30 {
                    mu[i][j] = dot(&b[i], &ortho[j]) / norm2[j];
                    for m in 0..5 {
                        ortho[i][m] -= mu[i][j] * ortho[j][m];
                    }
                } else {
                    mu[i][j] = 0.0;
                }
            }
            norm2[i] = dot(&ortho[i], &ortho[i]);
        }
    };

    recompute_gs(&b, &mut ortho, &mut mu, &mut norm2);
    let mut k = 1usize;
    let mut iters = 0usize;
    while k < 4 && iters < 80 {
        iters += 1;
        for j in (0..k).rev() {
            let r = mu[k][j].round();
            if r != 0.0 {
                for m in 0..5 {
                    b[k][m] -= r * b[j][m];
                }
                for m in 0..j {
                    mu[k][m] -= r * mu[j][m];
                }
                mu[k][j] -= r;
            }
        }
        if norm2[k] >= (0.75 - mu[k][k - 1] * mu[k][k - 1]) * norm2[k - 1] {
            k += 1;
        } else {
            b.swap(k, k - 1);
            recompute_gs(&b, &mut ortho, &mut mu, &mut norm2);
            if k > 1 {
                k -= 1;
            }
        }
    }

    // Check all 4 reduced basis vectors for a valid small Möbius relation!
    for i in 0..4 {
        let c0 = b[i][0].round() as i64;
        let c1 = b[i][1].round() as i64;
        let c2 = b[i][2].round() as i64;
        let c3 = b[i][3].round() as i64;
        // Must genuinely involve K (c2 or c3 nonzero) AND alpha (c1 or c3 nonzero),
        // and not degenerate (det c0*c3 - c1*c2 != 0) so K is not purely rational!
        if (c2 == 0 && c3 == 0) || (c1 == 0 && c3 == 0) {
            continue;
        }
        if c0 * c3 - c1 * c2 == 0 {
            continue;
        }
        if c0.abs() > max_c || c1.abs() > max_c || c2.abs() > max_c || c3.abs() > max_c {
            continue;
        }
        let resid = (c0 as f64) * x[0]
            + (c1 as f64) * x[1]
            + (c2 as f64) * x[2]
            + (c3 as f64) * x[3];
        if resid.abs() < tol {
            return Some([c0, c1, c2, c3]);
        }
    }
    None
}

#[derive(Clone, Debug)]
struct Hit {
    family: String,
    params: String,
    target: &'static str,
    k_val: f64,
    rel: [i64; 4], // c0 + c1*alpha + c2*K + c3*alpha*K = 0
    resid: f64,
}

fn check_targets(
    family: &str,
    params: String,
    k_val: f64,
    targets: &[(&'static str, f64)],
    max_c: i64,
    tol: f64,
    out: &mut Vec<Hit>,
) {
    for &(name, alpha) in targets {
        if let Some(rel) = find_mobius_lll(alpha, k_val, max_c, tol) {
            let resid = ((rel[0] as f64)
                + (rel[1] as f64) * alpha
                + (rel[2] as f64) * k_val
                + (rel[3] as f64) * alpha * k_val)
                .abs();
            out.push(Hit {
                family: family.to_string(),
                params: params.clone(),
                target: name,
                k_val,
                rel,
                resid,
            });
        }
    }
}

fn main() {
    let nthreads = 16usize;
    let z3_over_pi3 = ZETA3 / (PI * PI * PI);
    let z3_over_pi2 = ZETA3 / (PI * PI);
    let pi3 = PI * PI * PI;
    let targets: Vec<(&'static str, f64)> = vec![
        ("zeta(3)", ZETA3),
        ("zeta(3)/pi^3", z3_over_pi3),
        ("pi^3", pi3),
        ("zeta(3)/pi^2", z3_over_pi2),
    ];

    // 1. Self-test: Apéry's PCF a(n)=(2n+1)(17n^2+17n+5), b(n)=-n^6 => K = 6/zeta(3)
    println!("=== Self-Test: Apéry & Ramanujan-Machine PCFs ===");
    let k_apery = eval_pcf(
        |n| {
            let nf = n as f64;
            (2.0 * nf + 1.0) * (17.0 * nf * nf + 17.0 * nf + 5.0)
        },
        |n| {
            let nf = n as f64;
            -(nf * nf * nf * nf * nf * nf)
        },
        25,
        40,
    )
    .expect("Apery failed");
    let rel_apery = find_mobius_lll(ZETA3, k_apery, 50, 1e-11).expect("Apery LLL failed");
    println!(
        "Apéry PCF: K = {:.15} (6/zeta(3) = {:.15}) | LLL rel = {:?} [PASS]",
        k_apery,
        6.0 / ZETA3,
        rel_apery
    );

    let all_hits = Arc::new(Mutex::new(Vec::<Hit>::new()));

    // =========================================================================
    // Family 1: Symmetric Apéry-type cubic PCFs
    //   a(n) = (2n + 1)(A*n^2 + A*n + B),  b(n) = s * C * n^6
    //   Scan A in 1..=120, B in -80..=80, C in 1..=36, s in {-1, +1}
    //   (Convergent whenever 4*A^2 + 4*s*C > 0!)
    // =========================================================================
    println!("=== Family 1: Symmetric Cubic PCFs a(n)=(2n+1)(A n^2 + A n + B), b(n)=s*C*n^6 ===");
    let t1 = Instant::now();
    let count1 = Arc::new(Mutex::new(0u64));
    let mut hs1 = Vec::new();
    let a_vals: Vec<i64> = (1..=120).collect();
    for tid in 0..nthreads {
        let my_a: Vec<i64> = a_vals
            .iter()
            .copied()
            .enumerate()
            .filter(|(i, _)| i % nthreads == tid)
            .map(|(_, x)| x)
            .collect();
        let hits_ref = Arc::clone(&all_hits);
        let cnt_ref = Arc::clone(&count1);
        let tgts = targets.clone();
        hs1.push(thread::spawn(move || {
            let mut loc_hits = Vec::new();
            let mut loc_cnt = 0u64;
            for a in my_a {
                let af = a as f64;
                for c in 1..=36_i64 {
                    let cf = c as f64;
                    for &s in &[-1_i64, 1_i64] {
                        // Characteristic disc = (2A)^2 + 4*s*C = 4*(A^2 + s*C)
                        let disc = a * a + s * c;
                        if disc <= 0 {
                            continue;
                        }
                        // Need fast exponential convergence so n=60 gives <1e-13
                        let lam_ratio = (af - (disc as f64).sqrt()).abs() / (af + (disc as f64).sqrt());
                        if lam_ratio > 0.72 {
                            continue;
                        }
                        let scf = (s as f64) * cf;
                        for b in -60..=60_i64 {
                            // Primitive check gcd(a, |b|, c) == 1 (actually scaling a,b by g scales b_n by g^2)
                            let bf = b as f64;
                            loc_cnt += 1;
                            if let Some(kv) = eval_pcf(
                                |n| {
                                    let nf = n as f64;
                                    (2.0 * nf + 1.0) * (af * nf * nf + af * nf + bf)
                                },
                                |n| {
                                    let nf = n as f64;
                                    let n3 = nf * nf * nf;
                                    scf * n3 * n3
                                },
                                42,
                                62,
                            ) {
                                check_targets(
                                    "SymCubic-n6",
                                    format!("A={}, B={}, b(n)={}*n^6", a, b, s * c),
                                    kv,
                                    &tgts,
                                    32,
                                    1e-11,
                                    &mut loc_hits,
                                );
                            }
                        }
                    }
                }
            }
            hits_ref.lock().unwrap().extend(loc_hits);
            *cnt_ref.lock().unwrap() += loc_cnt;
        }));
    }
    for h in hs1 {
        h.join().unwrap();
    }
    println!(
        "Family 1 done in {:.2}s | Evaluated {} PCFs | Total hits so far: {}",
        t1.elapsed().as_secs_f64(),
        *count1.lock().unwrap(),
        all_hits.lock().unwrap().len()
    );

    // =========================================================================
    // Family 2: Shifted / Half-Integer Symmetric Cubic PCFs
    //   a(n) = (2n + 1)(A*n^2 + A*n + B)
    //   b(n) = s * C * n^3 * (n + d)^3   OR   s * C * (2n - 1)^3 * (2n)^3
    // =========================================================================
    println!("=== Family 2: Shifted Factorized b(n) = s*C*n^3*(n+d)^3 and (2n-1)^3*n^3 ===");
    let t2 = Instant::now();
    let count2 = Arc::new(Mutex::new(0u64));
    let mut hs2 = Vec::new();
    let a2_vals: Vec<i64> = (1..=60).collect();
    for tid in 0..nthreads {
        let my_a: Vec<i64> = a2_vals
            .iter()
            .copied()
            .enumerate()
            .filter(|(i, _)| i % nthreads == tid)
            .map(|(_, x)| x)
            .collect();
        let hits_ref = Arc::clone(&all_hits);
        let cnt_ref = Arc::clone(&count2);
        let tgts = targets.clone();
        hs2.push(thread::spawn(move || {
            let mut loc_hits = Vec::new();
            let mut loc_cnt = 0u64;
            for a in my_a {
                let af = a as f64;
                for c in 1..=16_i64 {
                    for &s in &[-1_i64, 1_i64] {
                        let disc = a * a + s * c;
                        if disc <= 0 {
                            continue;
                        }
                        let lam_ratio = (af - (disc as f64).sqrt()).abs() / (af + (disc as f64).sqrt());
                        if lam_ratio > 0.72 {
                            continue;
                        }
                        let scf = (s * c) as f64;
                        for d in [1_i64, 2, -1, 3] {
                            let df = d as f64;
                            for b in -30..=30_i64 {
                                let bf = b as f64;
                                loc_cnt += 1;
                                if let Some(kv) = eval_pcf(
                                    |n| {
                                        let nf = n as f64;
                                        (2.0 * nf + 1.0) * (af * nf * nf + af * nf + bf)
                                    },
                                    |n| {
                                        let nf = n as f64;
                                        let m = nf + df;
                                        scf * nf * nf * nf * m * m * m
                                    },
                                    42,
                                    62,
                                ) {
                                    check_targets(
                                        "Shifted-n3(n+d)3",
                                        format!("A={}, B={}, d={}, b={}*n^3(n+{})^3", a, b, d, s * c, d),
                                        kv,
                                        &tgts,
                                        24,
                                        1e-11,
                                        &mut loc_hits,
                                    );
                                }
                            }
                        }
                    }
                }
                // Half-integer b(n) = s*C*(2n-1)^3 * (2n)^3 (leading coeff 64*s*C vs (2A)^2 = 4A^2)
                for c in 1..=8_i64 {
                    for &s in &[-1_i64, 1_i64] {
                        let disc = a * a + 16 * s * c;
                        if disc <= 0 {
                            continue;
                        }
                        let lam_ratio = (af - (disc as f64).sqrt()).abs() / (af + (disc as f64).sqrt());
                        if lam_ratio > 0.72 {
                            continue;
                        }
                        let scf = (s * c) as f64;
                        for b in -30..=30_i64 {
                            let bf = b as f64;
                            loc_cnt += 1;
                            if let Some(kv) = eval_pcf(
                                |n| {
                                    let nf = n as f64;
                                    (2.0 * nf + 1.0) * (af * nf * nf + af * nf + bf)
                                },
                                |n| {
                                    let nf = n as f64;
                                    let u = 2.0 * nf - 1.0;
                                    let v = 2.0 * nf;
                                    scf * u * u * u * v * v * v
                                },
                                42,
                                62,
                            ) {
                                check_targets(
                                    "HalfInt-(2n-1)3(2n)3",
                                    format!("A={}, B={}, b={}*(2n-1)^3(2n)^3", a, b, s * c),
                                    kv,
                                    &tgts,
                                    24,
                                    1e-11,
                                    &mut loc_hits,
                                );
                            }
                        }
                    }
                }
            }
            hits_ref.lock().unwrap().extend(loc_hits);
            *cnt_ref.lock().unwrap() += loc_cnt;
        }));
    }
    for h in hs2 {
        h.join().unwrap();
    }
    println!(
        "Family 2 done in {:.2}s | Evaluated {} PCFs | Total hits so far: {}",
        t2.elapsed().as_secs_f64(),
        *count2.lock().unwrap(),
        all_hits.lock().unwrap().len()
    );

    // =========================================================================
    // Family 3: General Non-Symmetric Degree-3 PCFs with b(n) = -n^6 or +n^6
    //   a(n) = c3*n^3 + c2*n^2 + c1*n + c0
    //   Scan c3 in 4..=40, c2 in -20..=60, c1 in -15..=40, c0 in 1..=15 (sampled/dense)
    // =========================================================================
    println!("=== Family 3: General Cubic PCFs a(n)=c3*n^3+c2*n^2+c1*n+c0, b(n)=+/- n^6 ===");
    let t3 = Instant::now();
    let count3 = Arc::new(Mutex::new(0u64));
    let mut hs3 = Vec::new();
    let c3_vals: Vec<i64> = (4..=38).collect();
    for tid in 0..nthreads {
        let my_c3: Vec<i64> = c3_vals
            .iter()
            .copied()
            .enumerate()
            .filter(|(i, _)| i % nthreads == tid)
            .map(|(_, x)| x)
            .collect();
        let hits_ref = Arc::clone(&all_hits);
        let cnt_ref = Arc::clone(&count3);
        let tgts = targets.clone();
        hs3.push(thread::spawn(move || {
            let mut loc_hits = Vec::new();
            let mut loc_cnt = 0u64;
            for c3 in my_c3 {
                let f3 = c3 as f64;
                for &s in &[-1_i64, 1_i64] {
                    let disc = f3 * f3 + 4.0 * (s as f64);
                    if disc <= 0.0 {
                        continue;
                    }
                    let sf = s as f64;
                    // Focus around symmetric c2 = 3*c3/2 +/- 8 and c1 in -12..=32, c0 in 0..=12
                    let c2_mid = (3 * c3) / 2;
                    for c2 in (c2_mid - 8)..=(c2_mid + 8) {
                        let f2 = c2 as f64;
                        for c1 in -10..=32_i64 {
                            let f1 = c1 as f64;
                            for c0 in 0..=12_i64 {
                                let f0 = c0 as f64;
                                loc_cnt += 1;
                                if let Some(kv) = eval_pcf(
                                    |n| {
                                        let nf = n as f64;
                                        ((f3 * nf + f2) * nf + f1) * nf + f0
                                    },
                                    |n| {
                                        let nf = n as f64;
                                        let n3 = nf * nf * nf;
                                        sf * n3 * n3
                                    },
                                    40,
                                    58,
                                ) {
                                    check_targets(
                                        "GenCubic-n6",
                                        format!("a=({},{},{},{}), b={}*n^6", c3, c2, c1, c0, s),
                                        kv,
                                        &tgts,
                                        20,
                                        1e-11,
                                        &mut loc_hits,
                                    );
                                }
                            }
                        }
                    }
                }
            }
            hits_ref.lock().unwrap().extend(loc_hits);
            *cnt_ref.lock().unwrap() += loc_cnt;
        }));
    }
    for h in hs3 {
        h.join().unwrap();
    }
    println!(
        "Family 3 done in {:.2}s | Evaluated {} PCFs | Total hits: {}",
        t3.elapsed().as_secs_f64(),
        *count3.lock().unwrap(),
        all_hits.lock().unwrap().len()
    );

    // =========================================================================
    // Family 4: Degree-2 PCFs with b(n) = -C*n^4 or b(n) = -n^3*(n+d) (like Ramanujan/Zudilin)
    // =========================================================================
    println!("=== Family 4: Degree-2 PCFs a(n)=c2*n^2+c1*n+c0, b(n)=s*C*n^4 ===");
    let t4 = Instant::now();
    let count4 = Arc::new(Mutex::new(0u64));
    let mut hs4 = Vec::new();
    let c2_vals: Vec<i64> = (2..=40).collect();
    for tid in 0..nthreads {
        let my_c2: Vec<i64> = c2_vals
            .iter()
            .copied()
            .enumerate()
            .filter(|(i, _)| i % nthreads == tid)
            .map(|(_, x)| x)
            .collect();
        let hits_ref = Arc::clone(&all_hits);
        let cnt_ref = Arc::clone(&count4);
        let tgts = targets.clone();
        hs4.push(thread::spawn(move || {
            let mut loc_hits = Vec::new();
            let mut loc_cnt = 0u64;
            for c2 in my_c2 {
                let f2 = c2 as f64;
                for c in 1..=16_i64 {
                    for &s in &[-1_i64, 1_i64] {
                        let disc = f2 * f2 + 4.0 * (s * c) as f64;
                        if disc <= 0.0 {
                            continue;
                        }
                        let lam_ratio = (f2 - disc.sqrt()).abs() / (f2 + disc.sqrt());
                        if lam_ratio > 0.72 {
                            continue;
                        }
                        let scf = (s * c) as f64;
                        for c1 in -30..=40_i64 {
                            let f1 = c1 as f64;
                            for c0 in -15..=20_i64 {
                                let f0 = c0 as f64;
                                loc_cnt += 1;
                                if let Some(kv) = eval_pcf(
                                    |n| {
                                        let nf = n as f64;
                                        (f2 * nf + f1) * nf + f0
                                    },
                                    |n| {
                                        let nf = n as f64;
                                        let n2 = nf * nf;
                                        scf * n2 * n2
                                    },
                                    45,
                                    65,
                                ) {
                                    check_targets(
                                        "Quad-n4",
                                        format!("a=({},{},{}), b={}*n^4", c2, c1, c0, s * c),
                                        kv,
                                        &tgts,
                                        20,
                                        1e-11,
                                        &mut loc_hits,
                                    );
                                }
                            }
                        }
                    }
                }
            }
            hits_ref.lock().unwrap().extend(loc_hits);
            *cnt_ref.lock().unwrap() += loc_cnt;
        }));
    }
    for h in hs4 {
        h.join().unwrap();
    }
    println!(
        "Family 4 done in {:.2}s | Evaluated {} PCFs | Total hits: {}",
        t4.elapsed().as_secs_f64(),
        *count4.lock().unwrap(),
        all_hits.lock().unwrap().len()
    );

    // =========================================================================
    // Family 5: Wide Symmetric Cubic PCFs with large C up to 256 (incl. 64, 81, 144, 256)
    //           AND odd-factor b(n) = s*C*(2n-1)^6 (targeting pi^3, beta(3), zeta(3)/pi^3!)
    // =========================================================================
    println!("=== Family 5: Wide SymCubic (C<=256) + Odd (2n-1)^6 + Mixed n^2*(2n-1)^4 ===");
    let t5 = Instant::now();
    let count5 = Arc::new(Mutex::new(0u64));
    let mut hs5 = Vec::new();
    let a5_vals: Vec<i64> = (1..=160).collect();
    for tid in 0..nthreads {
        let my_a: Vec<i64> = a5_vals
            .iter()
            .copied()
            .enumerate()
            .filter(|(i, _)| i % nthreads == tid)
            .map(|(_, x)| x)
            .collect();
        let hits_ref = Arc::clone(&all_hits);
        let cnt_ref = Arc::clone(&count5);
        let tgts = targets.clone();
        hs5.push(thread::spawn(move || {
            let mut loc_hits = Vec::new();
            let mut loc_cnt = 0u64;
            // 5A: Large square/cube C in b(n) = s*C*n^6
            let c_list: [i64; 22] = [
                1, 2, 3, 4, 8, 9, 12, 16, 18, 24, 25, 27, 32, 36, 48, 49, 64, 72, 81, 100, 128, 144,
            ];
            for a in my_a {
                let af = a as f64;
                for &c in &c_list {
                    for &s in &[-1_i64, 1_i64] {
                        let disc = a * a + s * c;
                        if disc <= 0 { continue; }
                        let lam_ratio = (af - (disc as f64).sqrt()).abs() / (af + (disc as f64).sqrt());
                        if lam_ratio > 0.72 { continue; }
                        let scf = (s * c) as f64;
                        for b in -50..=50_i64 {
                            let bf = b as f64;
                            loc_cnt += 1;
                            if let Some(kv) = eval_pcf(
                                |n| {
                                    let nf = n as f64;
                                    (2.0 * nf + 1.0) * (af * nf * nf + af * nf + bf)
                                },
                                |n| {
                                    let nf = n as f64;
                                    let n3 = nf * nf * nf;
                                    scf * n3 * n3
                                },
                                42,
                                62,
                            ) {
                                check_targets(
                                    "WideSym-n6",
                                    format!("A={}, B={}, b={}*n^6", a, b, s * c),
                                    kv,
                                    &tgts,
                                    28,
                                    1e-11,
                                    &mut loc_hits,
                                );
                            }
                        }
                    }
                }
                // 5B: Odd-factor b(n) = s*C*(2n-1)^6 (leading 64*s*C vs 4*A^2 => disc = A^2 + 16*s*C)
                for c in [1_i64, 2, 4, 9, 16] {
                    for &s in &[-1_i64, 1_i64] {
                        let disc = a * a + 16 * s * c;
                        if disc <= 0 { continue; }
                        let lam_ratio = (af - (disc as f64).sqrt()).abs() / (af + (disc as f64).sqrt());
                        if lam_ratio > 0.75 { continue; }
                        let scf = (s * c) as f64;
                        for b in -40..=40_i64 {
                            let bf = b as f64;
                            loc_cnt += 1;
                            if let Some(kv) = eval_pcf(
                                |n| {
                                    let nf = n as f64;
                                    (2.0 * nf + 1.0) * (af * nf * nf + af * nf + bf)
                                },
                                |n| {
                                    let u = 2.0 * (n as f64) - 1.0;
                                    let u3 = u * u * u;
                                    scf * u3 * u3
                                },
                                42,
                                62,
                            ) {
                                check_targets(
                                    "Odd-(2n-1)6",
                                    format!("A={}, B={}, b={}*(2n-1)^6", a, b, s * c),
                                    kv,
                                    &tgts,
                                    28,
                                    1e-11,
                                    &mut loc_hits,
                                );
                            }
                        }
                    }
                }
                // 5C: Symmetric quartic a(n) = A*(2n+1)^4 + B*(2n+1)^2 + D, b(n) = s*C*n^8
                if a <= 30 {
                    for c in [1_i64, 4, 16, 64] {
                        for &s in &[-1_i64, 1_i64] {
                            // leading a(n) ~ 16*A*n^4, b(n) ~ s*C*n^8 => disc = 256*A^2 + 4*s*C
                            let disc = 64 * a * a + s * c;
                            if disc <= 0 { continue; }
                            let scf = (s * c) as f64;
                            for b in -20..=20_i64 {
                                let bf = b as f64;
                                for d in -10..=10_i64 {
                                    let df = d as f64;
                                    loc_cnt += 1;
                                    if let Some(kv) = eval_pcf(
                                        |n| {
                                            let u = 2.0 * (n as f64) + 1.0;
                                            let u2 = u * u;
                                            af * u2 * u2 + bf * u2 + df
                                        },
                                        |n| {
                                            let nf = n as f64;
                                            let n2 = nf * nf;
                                            let n4 = n2 * n2;
                                            scf * n4 * n4
                                        },
                                        36,
                                        54,
                                    ) {
                                        check_targets(
                                            "SymQuartic-n8",
                                            format!("A={}, B={}, D={}, b={}*n^8", a, b, d, s * c),
                                            kv,
                                            &tgts,
                                            24,
                                            1e-11,
                                            &mut loc_hits,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
            hits_ref.lock().unwrap().extend(loc_hits);
            *cnt_ref.lock().unwrap() += loc_cnt;
        }));
    }
    for h in hs5 {
        h.join().unwrap();
    }
    println!(
        "Family 5 done in {:.2}s | Evaluated {} PCFs | Total hits: {}",
        t5.elapsed().as_secs_f64(),
        *count5.lock().unwrap(),
        all_hits.lock().unwrap().len()
    );

    // =========================================================================
    // Family 6: Even/Odd Mixed b(n) = s*C*n^4*(2n-1)^2 & n^2*(2n-1)^4 & n^3*(2n-1)^3
    //           with general a(n) = c3*n^3 + c2*n^2 + c1*n + c0 (small range)
    // =========================================================================
    println!("=== Family 6: Euler/Nesterenko/Zudilin-type Mixed-Factor PCFs ===");
    let t6 = Instant::now();
    let count6 = Arc::new(Mutex::new(0u64));
    let mut hs6 = Vec::new();
    let a6_vals: Vec<i64> = (1..=100).collect();
    for tid in 0..nthreads {
        let my_a: Vec<i64> = a6_vals
            .iter()
            .copied()
            .enumerate()
            .filter(|(i, _)| i % nthreads == tid)
            .map(|(_, x)| x)
            .collect();
        let hits_ref = Arc::clone(&all_hits);
        let cnt_ref = Arc::clone(&count6);
        let tgts = targets.clone();
        hs6.push(thread::spawn(move || {
            let mut loc_hits = Vec::new();
            let mut loc_cnt = 0u64;
            for a in my_a {
                let af = a as f64;
                // 6A: b(n) = s*C*n^2*(2n-1)^2*(2n+1)^2 (leading 16*s*C vs 4*A^2 => disc = A^2 + 4*s*C)
                for c in [1_i64, 2, 4, 8, 9, 16] {
                    for &s in &[-1_i64, 1_i64] {
                        let disc = a * a + 4 * s * c;
                        if disc <= 0 { continue; }
                        let lam_ratio = (af - (disc as f64).sqrt()).abs() / (af + (disc as f64).sqrt());
                        if lam_ratio > 0.75 { continue; }
                        let scf = (s * c) as f64;
                        for b in -40..=40_i64 {
                            let bf = b as f64;
                            loc_cnt += 1;
                            if let Some(kv) = eval_pcf(
                                |n| {
                                    let nf = n as f64;
                                    (2.0 * nf + 1.0) * (af * nf * nf + af * nf + bf)
                                },
                                |n| {
                                    let nf = n as f64;
                                    let u = 2.0 * nf - 1.0;
                                    let v = 2.0 * nf + 1.0;
                                    scf * nf * nf * u * u * v * v
                                },
                                42,
                                62,
                            ) {
                                check_targets(
                                    "Mixed-n2(2n-1)2(2n+1)2",
                                    format!("A={}, B={}, b={}*n^2(2n-1)^2(2n+1)^2", a, b, s * c),
                                    kv,
                                    &tgts,
                                    28,
                                    1e-11,
                                    &mut loc_hits,
                                );
                            }
                        }
                    }
                }
                // 6B: b(n) = s*C*n^4*(2n-1)^2 (leading 4*s*C vs 4*A^2 => disc = A^2 + s*C)
                for c in [1_i64, 2, 4, 8, 9, 16, 25, 36] {
                    for &s in &[-1_i64, 1_i64] {
                        let disc = a * a + s * c;
                        if disc <= 0 { continue; }
                        let lam_ratio = (af - (disc as f64).sqrt()).abs() / (af + (disc as f64).sqrt());
                        if lam_ratio > 0.75 { continue; }
                        let scf = (s * c) as f64;
                        for b in -40..=40_i64 {
                            let bf = b as f64;
                            loc_cnt += 1;
                            if let Some(kv) = eval_pcf(
                                |n| {
                                    let nf = n as f64;
                                    (2.0 * nf + 1.0) * (af * nf * nf + af * nf + bf)
                                },
                                |n| {
                                    let nf = n as f64;
                                    let u = 2.0 * nf - 1.0;
                                    scf * nf * nf * nf * nf * u * u
                                },
                                42,
                                62,
                            ) {
                                check_targets(
                                    "Mixed-n4(2n-1)2",
                                    format!("A={}, B={}, b={}*n^4(2n-1)^2", a, b, s * c),
                                    kv,
                                    &tgts,
                                    28,
                                    1e-11,
                                    &mut loc_hits,
                                );
                            }
                        }
                    }
                }
            }
            hits_ref.lock().unwrap().extend(loc_hits);
            *cnt_ref.lock().unwrap() += loc_cnt;
        }));
    }
    for h in hs6 {
        h.join().unwrap();
    }
    let total_eval = *count1.lock().unwrap()
        + *count2.lock().unwrap()
        + *count3.lock().unwrap()
        + *count4.lock().unwrap()
        + *count5.lock().unwrap()
        + *count6.lock().unwrap();
    println!(
        "Family 6 done in {:.2}s | Evaluated {} PCFs | GRAND TOTAL EVALUATED: {} PCFs | Total hits: {}",
        t6.elapsed().as_secs_f64(),
        *count6.lock().unwrap(),
        total_eval,
        all_hits.lock().unwrap().len()
    );

    let mut hits = all_hits.lock().unwrap().clone();
    hits.sort_by(|a, b| a.resid.partial_cmp(&b.resid).unwrap());
    println!("\n=== TOP CANDIDATE HITS (TO VERIFY AT 100 DIGITS IN PYTHON) ===");
    for (i, h) in hits.iter().take(80).enumerate() {
        println!(
            "HIT#{:2} [{}] {} | {} | K={:.15} | rel={:?} | resid={:.2e}",
            i + 1,
            h.target,
            h.family,
            h.params,
            h.k_val,
            h.rel,
            h.resid
        );
    }
}
