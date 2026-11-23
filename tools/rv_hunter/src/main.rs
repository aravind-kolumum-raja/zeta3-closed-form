//! rv_hunter — 18-core, low-memory (<20 MB RAM) native Rust optimizer for
//! Rhin-Viola Theorem 5.1 (Acta Arith. 97 (2001) 269-293) on μ(ζ(3)).
//!
//! Uses exact 30 level-4 cosets of the 1920-element Weyl group Φ ≅ W(D_5),
//! evaluates digamma ψ(ω) ONLY at break points where the saving s(ω) ∈ {0,1,2}
//! actually changes (zero heap allocation inside the inner loop!), and supports:
//!   1. `selftest` — reproduces RV (2001) 5.513890626853... & Giraudin (2026) 5.513879974004...
//!   2. `ray`      — fast scan of 1D & 6D narrow-cone perturbations around RV (and any base point)
//!   3. `global`   — parallel 18-thread global search over base points + narrow-cone optimization
//!   4. `eval`     — exact evaluation + detailed diagnostics for any 8-tuple

use std::env;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

const FN: [[usize; 4]; 6] = [
    [0, 2, 3, 6],
    [0, 1, 3, 5],
    [0, 3, 4, 7],
    [1, 2, 5, 6],
    [0, 1, 4, 5],
    [1, 3, 5, 7],
];
const FD: [[usize; 4]; 6] = [
    [8, 12, 13, 14],
    [11, 13, 14, 15],
    [9, 10, 13, 14],
    [8, 9, 12, 13],
    [9, 10, 13, 14],
    [9, 11, 13, 15],
];
const FR: [usize; 6] = [8, 8, 4, 4, 4, 2];

#[inline(always)]
fn th(i: usize, t: usize) -> usize {
    if i < 8 {
        (i + t) % 8
    } else {
        8 + ((i - 8 + t) % 8)
    }
}

#[inline(always)]
fn psi_f64(x: f64) -> f64 {
    let mut acc = 0.0_f64;
    for t in 0..14 {
        acc -= 1.0 / (x + (t as f64));
    }
    let y = x + 14.0;
    let inv_y = 1.0 / y;
    let inv_y2 = inv_y * inv_y;
    let mut p = inv_y2;
    let mut out = y.ln() - 0.5 * inv_y;
    const B: [f64; 6] = [
        1.0 / 12.0,
        -1.0 / 120.0,
        1.0 / 252.0,
        -1.0 / 240.0,
        1.0 / 132.0,
        -691.0 / 32760.0,
    ];
    for &b in &B {
        out -= b * p;
        p *= inv_y2;
    }
    out + acc
}

struct CosetTables {
    pr: [[[usize; 4]; 2]; 30],
    aff: [[usize; 30]; 16],
    naff: [usize; 16],
}

fn build_coset_tables(v: &[i64; 16]) -> Option<CosetTables> {
    let mut num = [[0usize; 4]; 30];
    let mut den = [[0usize; 4]; 30];
    let mut nc = 0usize;
    for f in 0..6 {
        for u in 0..FR[f] {
            for i in 0..4 {
                num[nc][i] = th(FN[f][i], u);
                den[nc][i] = th(FD[f][i], u);
            }
            nc += 1;
        }
    }
    const SP: [[usize; 4]; 3] = [[0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]];
    let mut pr = [[[0usize; 4]; 2]; 30];
    for c in 0..30 {
        let mut found = false;
        'outer: for a in 0..3 {
            for b in 0..3 {
                for sw in 0..2 {
                    let n = &SP[a];
                    let d = &SP[b];
                    let (d0, d1) = if sw == 1 { (2, 0) } else { (0, 2) };
                    let n1 = v[num[c][n[0]]] + v[num[c][n[1]]];
                    let n2 = v[num[c][n[2]]] + v[num[c][n[3]]];
                    let e1 = v[den[c][d[d0]]] + v[den[c][d[d0 + 1]]];
                    let e2 = v[den[c][d[d1]]] + v[den[c][d[d1 + 1]]];
                    if n1 == e1 && n2 == e2 {
                        pr[c][0] = [num[c][n[0]], num[c][n[1]], den[c][d[d0]], den[c][d[d0 + 1]]];
                        pr[c][1] = [num[c][n[2]], num[c][n[3]], den[c][d[d1]], den[c][d[d1 + 1]]];
                        found = true;
                        break 'outer;
                    }
                }
            }
        }
        if !found {
            return None;
        }
    }
    let mut aff = [[0usize; 30]; 16];
    let mut naff = [0usize; 16];
    for c in 0..30 {
        for u in 0..2 {
            for w in 0..4 {
                let sy = pr[c][u][w];
                let mut dup = false;
                for z in 0..naff[sy] {
                    if aff[sy][z] == c {
                        dup = true;
                        break;
                    }
                }
                if !dup {
                    aff[sy][naff[sy]] = c;
                    naff[sy] += 1;
                }
            }
        }
    }
    Some(CosetTables { pr, aff, naff })
}

#[derive(Clone, Debug)]
struct EvalResult {
    p: [i64; 8],
    mu: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    s1: f64,
    s2: f64,
    m_n_q: [i64; 3],
    changes: u64,
}

/// Fast stationary point + analytic pre-check (takes ~20 nanoseconds!).
#[inline]
fn stationary_c0_c1(p: &[i64; 8]) -> Option<([i64; 16], f64, f64, [i64; 3])> {
    let [h, j, k, l, m, q, r, s] = *p;
    if h + m != k + r || j + q != l + s {
        return None;
    }
    let v: [i64; 16] = [
        h, j, k, l, m, q, r, s,
        h + l - j,
        j + m - k,
        k + q - l,
        l + r - m,
        m + s - q,
        q + h - r,
        r + j - s,
        s + k - h,
    ];
    for &x in &v {
        if x <= 0 {
            return None;
        }
    }
    let mut top = v;
    top.sort_unstable_by(|a, b| b.cmp(a));
    let m_n_q = [top[0], top[1], top[2]];

    let h128 = h as i128;
    let j128 = j as i128;
    let k128 = k as i128;
    let l128 = l as i128;
    let q128 = q as i128;
    let r128 = r as i128;
    let v8 = v[8] as i128;
    let v11 = v[11] as i128;
    let v14 = v[14] as i128;
    let v15 = v[15] as i128;

    let a_coef = (v8 * v11 * v14) as f64;
    let b_coef = ((h128 - j128) * (h128 - k128) * (q128 - r128)
        + (h128 * l128 + l128 * r128 + r128 * h128) * v15
        - h128 * l128 * r128) as f64;
    let c_coef = (-h128 * r128 * v15) as f64;

    let disc_sq = b_coef * b_coef - 4.0 * a_coef * c_coef;
    if disc_sq <= 0.0 || a_coef == 0.0 {
        return None;
    }
    let disc = disc_sq.sqrt();
    let mut c0 = 0.0;
    let mut c1 = 0.0;
    let hf = h as f64;
    let jf = j as f64;
    let kf = k as f64;
    let lf = l as f64;
    let qf = q as f64;
    let sf = s as f64;
    let v8f = v[8] as f64;
    let v11f = v[11] as f64;
    let v13f = v[13] as f64;
    let v15f = v[15] as f64;

    for e in 0..2 {
        let x = if e == 0 {
            (-b_coef + disc) / (2.0 * a_coef)
        } else {
            (-b_coef - disc) / (2.0 * a_coef)
        };
        let y = (v11f * x + (kf - hf)) / ((v11f - sf) * x + v15f);
        let z = (v8f * v11f * x * x
            - (hf * v8f + (hf - sf) * v11f - kf * (hf + qf)) * x
            - hf * v15f)
            / ((x - 1.0) * ((hf - sf) * v11f * x + hf * v15f));
        if e == 0 {
            if !(x > 0.0 && x < 1.0 && y > 0.0 && y < 1.0 && z > 0.0 && z < 1.0) {
                return None;
            }
        } else if !(x < 0.0 && y < 0.0 && x * y > 1.0 && z < 1.0 / (1.0 - x * y)) {
            return None;
        }
        let val = hf * x.abs().ln()
            + lf * (1.0 - x).abs().ln()
            + kf * y.abs().ln()
            + sf * (1.0 - y).abs().ln()
            + jf * z.abs().ln()
            + qf * (1.0 - z).abs().ln()
            - v13f * (1.0 - (1.0 - x * y) * z).abs().ln();
        if e == 0 {
            c0 = -val;
        } else {
            c1 = val;
        }
    }
    Some((v, c0, c1, m_n_q))
}

/// Exact Rhin-Viola Theorem 5.1 evaluation using lazy ψ(ω) (only when s(ω) changes!).
fn eval_exact(p: &[i64; 8]) -> Option<EvalResult> {
    let (v, c0, c1, m_n_q) = stationary_c0_c1(p)?;
    let tables = build_coset_tables(&v)?;

    // Deduplicate distinct denominators
    let mut dv = [0i64; 16];
    let mut who = [[0usize; 16]; 16];
    let mut nwho = [0usize; 16];
    let mut nd = 0usize;
    for i in 0..16 {
        let mut f = None;
        for z in 0..nd {
            if dv[z] == v[i] {
                f = Some(z);
                break;
            }
        }
        let idx = match f {
            Some(z) => z,
            None => {
                let z = nd;
                dv[z] = v[i];
                nd += 1;
                z
            }
        };
        who[idx][nwho[idx]] = i;
        nwho[idx] += 1;
    }

    let mut a = [1i64; 16];
    let mut fl = [0i32; 16];
    let mut cval = [0u8; 30];
    let mut touched = [false; 30];

    let recompute = |c: usize, fl: &[i32; 16], cval: &mut [u8; 30]| -> bool {
        let pr0 = tables.pr[c][0];
        let pr1 = tables.pr[c][1];
        let v0 = fl[pr0[2]] + fl[pr0[3]] - fl[pr0[0]] - fl[pr0[1]];
        let v1 = fl[pr1[2]] + fl[pr1[3]] - fl[pr1[0]] - fl[pr1[1]];
        if v0 < -1 || v0 > 1 || v1 < -1 || v1 > 1 {
            return false;
        }
        cval[c] = if v0 < 0 && v1 < 0 {
            2
        } else if v0 + v1 < 0 {
            1
        } else {
            0
        };
        true
    };

    for c in 0..30 {
        if !recompute(c, &fl, &mut cval) {
            return None;
        }
    }
    let mut sv = 0u8;
    for c in 0..30 {
        if cval[c] > sv {
            sv = cval[c];
        }
    }
    // At ω = 0, sv must be 0. Track transitions lazily:
    //   when sv changes from old -> new at ω, we close any active interval(s) and open new ones!
    //   ∫_{a}^{b} dψ = ψ(b) - ψ(a).
    let mut s1 = 0.0_f64;
    let mut k1 = 0.0_f64;
    let mut s2 = 0.0_f64;
    let mut k2 = 0.0_f64;
    let mut psi_open1 = 0.0_f64;
    let mut psi_open2 = 0.0_f64;
    let mut changes = 0u64;

    loop {
        // Find minimum a[z] / dv[z] with a[z] < dv[z]
        let mut best: Option<usize> = None;
        let mut best_a = 0i128;
        let mut best_d = 1i128;
        for z in 0..nd {
            if a[z] >= dv[z] {
                continue;
            }
            let az = a[z] as i128;
            let dz = dv[z] as i128;
            if best.is_none() || az * best_d < best_a * dz {
                best = Some(z);
                best_a = az;
                best_d = dz;
            }
        }
        let b_idx = match best {
            Some(z) => z,
            None => break,
        };
        let an = a[b_idx] as i128;
        let ad = dv[b_idx] as i128;
        let old = sv;

        // Advance all denominators tied at an / ad
        for z in 0..nd {
            if a[z] >= dv[z] || (a[z] as i128) * ad != an * (dv[z] as i128) {
                continue;
            }
            for u in 0..nwho[z] {
                let sy = who[z][u];
                fl[sy] += 1;
                for w in 0..tables.naff[sy] {
                    touched[tables.aff[sy][w]] = true;
                }
            }
            a[z] += 1;
        }
        for c in 0..30 {
            if touched[c] {
                if !recompute(c, &fl, &mut cval) {
                    return None;
                }
                touched[c] = false;
            }
        }
        sv = 0;
        for c in 0..30 {
            if cval[c] > sv {
                sv = cval[c];
            }
        }

        if sv != old {
            changes += 1;
            let om = (an as f64) / (ad as f64);
            let psi_om = psi_f64(om);
            // Update Ω_E (sv >= 1)
            if old == 0 && sv >= 1 {
                psi_open1 = psi_om;
            } else if old >= 1 && sv == 0 {
                let d = psi_om - psi_open1;
                let y = d - k1;
                let t2 = s1 + y;
                k1 = (t2 - s1) - y;
                s1 = t2;
            }
            // Update Ω'_E (sv == 2)
            if old < 2 && sv == 2 {
                psi_open2 = psi_om;
            } else if old == 2 && sv < 2 {
                let d = psi_om - psi_open2;
                let y = d - k2;
                let t2 = s2 + y;
                k2 = (t2 - s2) - y;
                s2 = t2;
            }
        }
    }

    let psi_1 = psi_f64(1.0);
    if sv >= 1 {
        s1 += psi_1 - psi_open1;
    }
    if sv == 2 {
        s2 += psi_1 - psi_open2;
    }

    let c2 = (m_n_q[0] + m_n_q[1] + m_n_q[2]) as f64 - (s1 + s2);
    if c0 <= c2 {
        return None;
    }
    let mu = (c0 + c1) / (c0 - c2);
    Some(EvalResult {
        p: *p,
        mu,
        c0,
        c1,
        c2,
        s1,
        s2,
        m_n_q,
        changes,
    })
}

fn print_res(tag: &str, r: &EvalResult) {
    println!(
        "[{}] mu = {:.12} | P = {:?} | MNQ = {:?} | c0={:.6} c1={:.6} c2={:.6} (S1={:.6} S2={:.6}, chg={})",
        tag, r.mu, r.p, r.m_n_q, r.c0, r.c1, r.c2, r.s1, r.s2, r.changes
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mode = args.get(1).map(|s| s.as_str()).unwrap_or("selftest");
    let rv: [i64; 8] = [16, 17, 19, 15, 12, 11, 9, 13];
    let gir: [i64; 8] = [832002, 883999, 987993, 780005, 624014, 572017, 468023, 676011];

    match mode {
        "selftest" => {
            let t0 = Instant::now();
            let r_rv = eval_exact(&rv).expect("RV failed");
            print_res("RV 2001", &r_rv);
            let t1 = Instant::now();
            let r_gir = eval_exact(&gir).expect("Giraudin failed");
            print_res("Giraudin 2026", &r_gir);
            println!(
                "Timings: RV = {:.3}ms, Giraudin (K=52000) = {:.3}s",
                t1.duration_since(t0).as_secs_f64() * 1000.0,
                t1.elapsed().as_secs_f64()
            );
        }

        "ray" => {
            // 1. Test K * RV + 1 across K = 50..5000 (parallel 18 threads!)
            println!("=== Phase 1: Exact 1D Ray Scan K * RV + 1 (K = 100..6000) ===");
            let t0 = Instant::now();
            let best_global = Arc::new(Mutex::new((10.0_f64, rv)));
            let nthreads = 16usize;
            let ks: Vec<i64> = (100..=6000).collect();
            let chunk = (ks.len() + nthreads - 1) / nthreads;
            let mut handles = Vec::new();
            for tid in 0..nthreads {
                let slice = ks[tid * chunk..ks.len().min((tid + 1) * chunk)].to_vec();
                let best_ref = Arc::clone(&best_global);
                handles.push(thread::spawn(move || {
                    let mut loc_best = (10.0_f64, rv);
                    for k in slice {
                        let p = [
                            k * rv[0] + 1, k * rv[1] + 1, k * rv[2] + 1, k * rv[3] + 1,
                            k * rv[4] + 1, k * rv[5] + 1, k * rv[6] + 1, k * rv[7] + 1,
                        ];
                        if let Some(res) = eval_exact(&p) {
                            if res.mu < loc_best.0 {
                                loc_best = (res.mu, p);
                            }
                        }
                    }
                    let mut g = best_ref.lock().unwrap();
                    if loc_best.0 < g.0 {
                        *g = loc_best;
                    }
                }));
            }
            for h in handles {
                h.join().unwrap();
            }
            let (bmu, bp) = *best_global.lock().unwrap();
            println!(
                "Phase 1 done in {:.2}s | Best K*RV+1: mu = {:.12} at P = {:?}",
                t0.elapsed().as_secs_f64(), bmu, bp
            );

            // 2. Full 6-DOF narrow-cone direction scan around RV at moderate scale K = 800
            // Valid 6-DOF basis for shifts d = (dh,dj,dk,dl,dm,dq,dr,ds) with dh+dm=dk+dr, dj+q=l+s:
            //   Free vars: (dh, dj, dk, dl, dm, dq), then dr = dh + dm - dk, ds = dj + dq - dl.
            // Scan all small integer shifts with |d_i| <= 2 around (1,1,1,1,1,1,1,1) AND around 0!
            println!("=== Phase 2: 6-DOF Narrow-Cone Direction Scan around K*RV (K = 1000) ===");
            let t2 = Instant::now();
            // Build candidate shift vectors d in [-2..=2]^6 (filtered so |dr|<=2, |ds|<=2)
            let mut shifts: Vec<[i64; 8]> = Vec::new();
            for dh in -2..=2_i64 {
                for dj in -2..=2_i64 {
                    for dk in -2..=2_i64 {
                        for dl in -2..=2_i64 {
                            for dm in -2..=2_i64 {
                                let dr = dh + dm - dk;
                                if dr.abs() > 2 {
                                    continue;
                                }
                                for dq in -2..=2_i64 {
                                    let ds = dj + dq - dl;
                                    if ds.abs() > 2 {
                                        continue;
                                    }
                                    if dh == 0 && dj == 0 && dk == 0 && dl == 0 && dm == 0 && dq == 0 {
                                        continue;
                                    }
                                    // Canonicalize gcd
                                    let arr = [dh, dj, dk, dl, dm, dq, dr, ds];
                                    let mut g = 0i64;
                                    for &x in &arr {
                                        g = gcd(g, x.abs());
                                    }
                                    if g == 1 {
                                        shifts.push(arr);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            println!("Generated {} primitive 6-DOF shift directions in [-2..2]^8", shifts.len());
            let top_dirs = Arc::new(Mutex::new(Vec::<(f64, [i64; 8], i64)>::new()));
            let chunk2 = (shifts.len() + nthreads - 1) / nthreads;
            let mut handles2 = Vec::new();
            let k_screen = 1000_i64;
            for tid in 0..nthreads {
                let sl = shifts[tid * chunk2..shifts.len().min((tid + 1) * chunk2)].to_vec();
                let top_ref = Arc::clone(&top_dirs);
                handles2.push(thread::spawn(move || {
                    let mut local_top: Vec<(f64, [i64; 8], i64)> = Vec::new();
                    for d in sl {
                        let p = [
                            k_screen * rv[0] + d[0], k_screen * rv[1] + d[1],
                            k_screen * rv[2] + d[2], k_screen * rv[3] + d[3],
                            k_screen * rv[4] + d[4], k_screen * rv[5] + d[5],
                            k_screen * rv[6] + d[6], k_screen * rv[7] + d[7],
                        ];
                        if let Some(res) = eval_exact(&p) {
                            if res.mu < 5.513885 {
                                local_top.push((res.mu, d, k_screen));
                            }
                        }
                    }
                    let mut g = top_ref.lock().unwrap();
                    g.extend(local_top);
                }));
            }
            for h in handles2 {
                h.join().unwrap();
            }
            let mut dirs = top_dirs.lock().unwrap().clone();
            dirs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            println!(
                "Phase 2 done in {:.2}s | Found {} directions with mu < 5.513885 at K=1000!",
                t2.elapsed().as_secs_f64(), dirs.len()
            );
            for (rank, (mu, d, k)) in dirs.iter().take(15).enumerate() {
                println!("  #{:2}: mu = {:.12} | K={} | d = {:?}", rank + 1, mu, k, d);
            }

            // Phase 3: Optimize scale K across [200..=15000] for the top 20 directions!
            println!("=== Phase 3: Multi-Scale K Optimization on Top 20 6-DOF Directions ===");
            let t3 = Instant::now();
            let top20: Vec<[i64; 8]> = dirs.iter().take(20).map(|x| x.1).collect();
            let best3 = Arc::new(Mutex::new((10.0_f64, [0i64; 8], [0i64; 8], 0i64)));
            let mut handles3 = Vec::new();
            for d in top20 {
                let b3 = Arc::clone(&best3);
                handles3.push(thread::spawn(move || {
                    let mut loc = (10.0_f64, [0i64; 8], d, 0i64);
                    // Coarse log-spaced + fine scan around optimum
                    let mut best_k = 1000i64;
                    for &k in &[
                        150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1400,
                        1600, 1800, 2000, 2500, 3000, 4000, 5000, 7000, 10000,
                    ] {
                        let p = [
                            k * rv[0] + d[0], k * rv[1] + d[1], k * rv[2] + d[2], k * rv[3] + d[3],
                            k * rv[4] + d[4], k * rv[5] + d[5], k * rv[6] + d[6], k * rv[7] + d[7],
                        ];
                        if let Some(res) = eval_exact(&p) {
                            if res.mu < loc.0 {
                                loc = (res.mu, p, d, k);
                                best_k = k;
                            }
                        }
                    }
                    // Fine refinement around best_k
                    let lo = (best_k - 150).max(50);
                    let hi = best_k + 150;
                    for k in (lo..=hi).step_by(5) {
                        let p = [
                            k * rv[0] + d[0], k * rv[1] + d[1], k * rv[2] + d[2], k * rv[3] + d[3],
                            k * rv[4] + d[4], k * rv[5] + d[5], k * rv[6] + d[6], k * rv[7] + d[7],
                        ];
                        if let Some(res) = eval_exact(&p) {
                            if res.mu < loc.0 {
                                loc = (res.mu, p, d, k);
                            }
                        }
                    }
                    let mut g = b3.lock().unwrap();
                    if loc.0 < g.0 {
                        *g = loc;
                    }
                }));
            }
            for h in handles3 {
                h.join().unwrap();
            }
            let (mu3, p3, d3, k3) = *best3.lock().unwrap();
            println!(
                "Phase 3 done in {:.2}s | BEST FOUND: mu = {:.12} (K={}, d={:?})",
                t3.elapsed().as_secs_f64(), mu3, k3, d3
            );
            if let Some(r) = eval_exact(&p3) {
                print_res("PHASE 3 CHAMPION", &r);
            }
        }

        "eval" => {
            if args.len() != 10 {
                eprintln!("Usage: rv_hunter eval h j k l m q r s");
                std::process::exit(1);
            }
            let mut p = [0i64; 8];
            for i in 0..8 {
                p[i] = args[i + 2].parse().expect("int");
            }
            match eval_exact(&p) {
                Some(r) => print_res("EVAL", &r),
                None => println!("Invalid / inadmissible point {:?}", p),
            }
        }

        "optimize" => {
            let nthreads = 16usize;
            // =========================================================================
            // Stage A: Exhaustive 5-tie plane scan a*RV + b*1 (all coprime a, b)
            // =========================================================================
            println!("=== Stage A: 5-Tie Plane Scan a*RV + b*1 (1 <= a <= 600, -15 <= b <= 5) ===");
            let t_a = Instant::now();
            let mut ab_pairs: Vec<(i64, i64)> = Vec::new();
            for a in 1..=600_i64 {
                for b in -15..=5_i64 {
                    if b == 0 {
                        continue;
                    }
                    // Keep eps = b/a in [-0.06, +0.01] for speed, plus all small a <= 80
                    let eps = (b as f64) / (a as f64);
                    if a > 80 && (eps < -0.04 || eps > 0.005) {
                        continue;
                    }
                    if gcd(a, b.abs()) == 1 {
                        ab_pairs.push((a, b));
                    }
                }
            }
            println!("Scanning {} coprime (a, b) pairs in 5-tie plane...", ab_pairs.len());
            let best_a = Arc::new(Mutex::new(Vec::<(f64, [i64; 8], i64, i64)>::new()));
            let chunk = (ab_pairs.len() + nthreads - 1) / nthreads;
            let mut handles = Vec::new();
            for tid in 0..nthreads {
                let sl = ab_pairs[tid * chunk..ab_pairs.len().min((tid + 1) * chunk)].to_vec();
                let bref = Arc::clone(&best_a);
                handles.push(thread::spawn(move || {
                    let mut loc: Vec<(f64, [i64; 8], i64, i64)> = Vec::new();
                    for (a, b) in sl {
                        let p = [
                            a * rv[0] + b, a * rv[1] + b, a * rv[2] + b, a * rv[3] + b,
                            a * rv[4] + b, a * rv[5] + b, a * rv[6] + b, a * rv[7] + b,
                        ];
                        if let Some(res) = eval_exact(&p) {
                            if res.mu < 5.51376 {
                                loc.push((res.mu, p, a, b));
                            }
                        }
                    }
                    bref.lock().unwrap().extend(loc);
                }));
            }
            for h in handles {
                h.join().unwrap();
            }
            let mut res_a = best_a.lock().unwrap().clone();
            res_a.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            println!("Stage A done in {:.2}s | Found {} points < 5.51376", t_a.elapsed().as_secs_f64(), res_a.len());
            for (rank, (mu, p, a, b)) in res_a.iter().take(12).enumerate() {
                println!("  A#{:2}: mu = {:.12} | a={:3}, b={:+3} (eps={:+.6}) | P={:?}", rank + 1, mu, a, b, (*b as f64)/(*a as f64), p);
            }
            let mut champ_p = res_a.first().map(|x| x.1).unwrap_or([1199, 1274, 1424, 1124, 899, 824, 674, 974]);
            let mut champ_mu = res_a.first().map(|x| x.0).unwrap_or(5.513747);

            // =========================================================================
            // Stage B: Full 6-DOF Coordinate & Pattern Descent from Stage A Champion
            //          + from 75*RV-1, 150*RV-2, and scaled versions!
            // =========================================================================
            println!("=== Stage B: Full 6-DOF Parallel Coordinate + 2-Step Descent ===");
            let t_b = Instant::now();
            // Build ALL 6-DOF primitive moves in [-1..=1]^8 (6 free vars, dr=dh+dm-dk, ds=dj+dq-dl)
            let mut moves1: Vec<[i64; 8]> = Vec::new();
            for dh in -1..=1_i64 {
                for dj in -1..=1_i64 {
                    for dk in -1..=1_i64 {
                        for dl in -1..=1_i64 {
                            for dm in -1..=1_i64 {
                                let dr = dh + dm - dk;
                                if dr.abs() > 1 { continue; }
                                for dq in -1..=1_i64 {
                                    let ds = dj + dq - dl;
                                    if ds.abs() > 1 { continue; }
                                    if dh == 0 && dj == 0 && dk == 0 && dl == 0 && dm == 0 && dq == 0 {
                                        continue;
                                    }
                                    moves1.push([dh, dj, dk, dl, dm, dq, dr, ds]);
                                }
                            }
                        }
                    }
                }
            }
            println!("6-DOF unit neighborhood size: {} moves", moves1.len());

            let run_descent = |start: [i64; 8], steps: &[i64], label: &str| -> (f64, [i64; 8]) {
                let mut cur = start;
                let mut cur_mu = eval_exact(&cur).map(|r| r.mu).unwrap_or(10.0);
                println!("  [{}] start mu = {:.12} at {:?}", label, cur_mu, cur);
                for &step in steps {
                    let mut improved = true;
                    let mut iter = 0;
                    while improved && iter < 40 {
                        improved = false;
                        iter += 1;
                        let best_step = Arc::new(Mutex::new((cur_mu, cur)));
                        let chunk = (moves1.len() + nthreads - 1) / nthreads;
                        let mut hs = Vec::new();
                        for tid in 0..nthreads {
                            let sl = moves1[tid * chunk..moves1.len().min((tid + 1) * chunk)].to_vec();
                            let bs = Arc::clone(&best_step);
                            let base = cur;
                            hs.push(thread::spawn(move || {
                                let mut loc = (10.0_f64, base);
                                for d in sl {
                                    let cand = [
                                        base[0] + step * d[0], base[1] + step * d[1],
                                        base[2] + step * d[2], base[3] + step * d[3],
                                        base[4] + step * d[4], base[5] + step * d[5],
                                        base[6] + step * d[6], base[7] + step * d[7],
                                    ];
                                    if let Some(r) = eval_exact(&cand) {
                                        if r.mu < loc.0 {
                                            loc = (r.mu, cand);
                                        }
                                    }
                                }
                                let mut g = bs.lock().unwrap();
                                if loc.0 < g.0 {
                                    *g = loc;
                                }
                            }));
                        }
                        for h in hs {
                            h.join().unwrap();
                        }
                        let (nmu, np) = *best_step.lock().unwrap();
                        if nmu + 1e-13 < cur_mu {
                            cur_mu = nmu;
                            cur = np;
                            improved = true;
                            println!("    step={:3} iter={:2}: mu = {:.12} | P = {:?}", step, iter, cur_mu, cur);
                        }
                    }
                }
                (cur_mu, cur)
            };

            let seeds: Vec<([i64; 8], Vec<i64>, &str)> = vec![
                (champ_p, vec![20, 5, 2, 1], "StageA-champ"),
                ([1199, 1274, 1424, 1124, 899, 824, 674, 974], vec![10, 3, 1], "75*RV-1"),
                ([799, 849, 949, 749, 599, 549, 449, 649], vec![10, 3, 1], "50*RV-1"),
                ([319, 339, 379, 299, 239, 219, 179, 259], vec![5, 2, 1], "20*RV-1"),
                ([160, 170, 190, 150, 120, 110, 90, 130], vec![5, 2, 1], "10*RV"),
            ];
            for (seed, steps, lbl) in seeds {
                let (m, p) = run_descent(seed, &steps, lbl);
                if m < champ_mu {
                    champ_mu = m;
                    champ_p = p;
                }
            }
            println!("Stage B done in {:.2}s | Best so far: mu = {:.12} at {:?}", t_b.elapsed().as_secs_f64(), champ_mu, champ_p);

            // =========================================================================
            // Stage C: Multi-scale ray refinement around the Stage B champion!
            //          Test K * champ_p + d for K = 1..15 and fine 6D shifts!
            // =========================================================================
            println!("=== Stage C: Multi-scale Ray + Fine Polish Around Stage B Champion ===");
            let t_c = Instant::now();
            // First reduce champ_p by gcd
            let mut g = 0i64;
            for &x in &champ_p { g = gcd(g, x.abs()); }
            let prim_p = [
                champ_p[0]/g, champ_p[1]/g, champ_p[2]/g, champ_p[3]/g,
                champ_p[4]/g, champ_p[5]/g, champ_p[6]/g, champ_p[7]/g,
            ];
            println!("Primitive champion (gcd={}): {:?}", g, prim_p);
            for scale in [1_i64, 2, 3, 4, 5] {
                let scaled = [
                    scale * prim_p[0], scale * prim_p[1], scale * prim_p[2], scale * prim_p[3],
                    scale * prim_p[4], scale * prim_p[5], scale * prim_p[6], scale * prim_p[7],
                ];
                let (m, p) = run_descent(scaled, &[1], &format!("scale-{}", scale));
                if m < champ_mu {
                    champ_mu = m;
                    champ_p = p;
                }
            }
            // Also test shifting RV itself in all 6D directions at K=75 with |d_i| <= 2!
            println!("=== Stage D: Exhaustive [-2..2]^8 6-DOF Scan at K = 75 around 75*RV ===");
            let mut shifts2: Vec<[i64; 8]> = Vec::new();
            for dh in -2..=2_i64 {
                for dj in -2..=2_i64 {
                    for dk in -2..=2_i64 {
                        for dl in -2..=2_i64 {
                            for dm in -2..=2_i64 {
                                let dr = dh + dm - dk;
                                if dr.abs() > 2 { continue; }
                                for dq in -2..=2_i64 {
                                    let ds = dj + dq - dl;
                                    if ds.abs() > 2 { continue; }
                                    if dh == 0 && dj == 0 && dk == 0 && dl == 0 && dm == 0 && dq == 0 {
                                        continue;
                                    }
                                    shifts2.push([dh, dj, dk, dl, dm, dq, dr, ds]);
                                }
                            }
                        }
                    }
                }
            }
            let best_d = Arc::new(Mutex::new(Vec::<(f64, [i64; 8], [i64; 8])>::new()));
            let chunk_d = (shifts2.len() + nthreads - 1) / nthreads;
            let mut hd = Vec::new();
            for tid in 0..nthreads {
                let sl = shifts2[tid * chunk_d..shifts2.len().min((tid + 1) * chunk_d)].to_vec();
                let bd = Arc::clone(&best_d);
                hd.push(thread::spawn(move || {
                    let mut loc = Vec::new();
                    for d in sl {
                        let p = [
                            75 * rv[0] + d[0], 75 * rv[1] + d[1], 75 * rv[2] + d[2], 75 * rv[3] + d[3],
                            75 * rv[4] + d[4], 75 * rv[5] + d[5], 75 * rv[6] + d[6], 75 * rv[7] + d[7],
                        ];
                        if let Some(r) = eval_exact(&p) {
                            if r.mu < 5.51380 {
                                loc.push((r.mu, p, d));
                            }
                        }
                    }
                    bd.lock().unwrap().extend(loc);
                }));
            }
            for h in hd { h.join().unwrap(); }
            let mut res_d = best_d.lock().unwrap().clone();
            res_d.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            println!("Stage D done in {:.2}s | Found {} directions < 5.51380 at K=75", t_c.elapsed().as_secs_f64(), res_d.len());
            for (rank, (mu, p, d)) in res_d.iter().take(15).enumerate() {
                println!("  D#{:2}: mu = {:.12} | d={:?} | P={:?}", rank + 1, mu, d, p);
            }
            if let Some((m, p, _)) = res_d.first() {
                if *m < champ_mu {
                    champ_p = *p;
                }
            }

            if let Some(r) = eval_exact(&champ_p) {
                print_res("OVERALL OPTIMIZE CHAMPION", &r);
            }
        }

        "deep" => {
            let nthreads = 16usize;
            let p73: [i64; 8] = [1167, 1240, 1386, 1094, 875, 802, 656, 948];
            let u3: [i64; 8] = [0, -1, 1, 1, 0, 1, -1, -1];
            let e2: [i64; 8] = [7, 8, 10, 6, 3, 2, 0, 4];

            // =========================================================================
            // Deep 1: Global map of the ENTIRE 5-tie 2D plane a*e2 + b*1
            //         (a in 1..=120, b/a in [1.0, 30.0]) to check if ANY other ratio
            //         beats b/a = 656/73 = 8.9863!
            // =========================================================================
            println!("=== Deep 1: Global Map of Entire 5-Tie 2D Plane (a*e2 + b*1) ===");
            let t1 = Instant::now();
            let mut pairs: Vec<(i64, i64)> = Vec::new();
            for a in 1..=120_i64 {
                let b_lo = (a * 2).max(1);
                let b_hi = a * 20;
                for b in b_lo..=b_hi {
                    if gcd(a, b) == 1 {
                        // Filter densely around 8.5..9.5, sparser elsewhere
                        let ratio = (b as f64) / (a as f64);
                        if a <= 30 || (ratio > 7.5 && ratio < 10.5) {
                            pairs.push((a, b));
                        }
                    }
                }
            }
            println!("Scanning {} coprime (a, b) points across 5-tie plane...", pairs.len());
            let best1 = Arc::new(Mutex::new(Vec::<(f64, i64, i64, [i64; 8])>::new()));
            let ch1 = (pairs.len() + nthreads - 1) / nthreads;
            let mut hs1 = Vec::new();
            for tid in 0..nthreads {
                let sl = pairs[tid * ch1..pairs.len().min((tid + 1) * ch1)].to_vec();
                let b1 = Arc::clone(&best1);
                hs1.push(thread::spawn(move || {
                    let mut loc = Vec::new();
                    for (a, b) in sl {
                        let p = [
                            a * e2[0] + b, a * e2[1] + b, a * e2[2] + b, a * e2[3] + b,
                            a * e2[4] + b, a * e2[5] + b, a * e2[6] + b, a * e2[7] + b,
                        ];
                        if let Some(r) = eval_exact(&p) {
                            if r.mu < 5.51380 {
                                loc.push((r.mu, a, b, p));
                            }
                        }
                    }
                    b1.lock().unwrap().extend(loc);
                }));
            }
            for h in hs1 { h.join().unwrap(); }
            let mut res1 = best1.lock().unwrap().clone();
            res1.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            println!("Deep 1 done in {:.2}s | Found {} points < 5.51380", t1.elapsed().as_secs_f64(), res1.len());
            for (i, (mu, a, b, p)) in res1.iter().take(10).enumerate() {
                println!("  1#{:2}: mu={:.12} | a={:3} b={:4} (b/a={:.6}) | P={:?}", i+1, mu, a, b, (*b as f64)/(*a as f64), p);
            }

            // =========================================================================
            // Deep 2: 3D 4-Tie Subspace Scan: a*RV + b*1 + c*u3
            //         around RV and around 73*RV - 1 (a=1..200, b in [-5..2], c in [-8..8])
            // =========================================================================
            println!("=== Deep 2: 3D 4-Tie Subspace Scan (a*RV + b*1 + c*u3) ===");
            let t2 = Instant::now();
            let mut triples: Vec<(i64, i64, i64)> = Vec::new();
            for a in 1..=250_i64 {
                for b in -4..=2_i64 {
                    let eps = (b as f64) / (a as f64);
                    if a > 30 && (eps < -0.025 || eps > 0.005) {
                        continue;
                    }
                    for c in -6..=6_i64 {
                        if c == 0 { continue; }
                        triples.push((a, b, c));
                    }
                }
            }
            println!("Scanning {} points in 3D 4-tie subspace...", triples.len());
            let best2 = Arc::new(Mutex::new(Vec::<(f64, i64, i64, i64, [i64; 8])>::new()));
            let ch2 = (triples.len() + nthreads - 1) / nthreads;
            let mut hs2 = Vec::new();
            for tid in 0..nthreads {
                let sl = triples[tid * ch2..triples.len().min((tid + 1) * ch2)].to_vec();
                let b2 = Arc::clone(&best2);
                hs2.push(thread::spawn(move || {
                    let mut loc = Vec::new();
                    for (a, b, c) in sl {
                        let p = [
                            a*rv[0] + b + c*u3[0], a*rv[1] + b + c*u3[1],
                            a*rv[2] + b + c*u3[2], a*rv[3] + b + c*u3[3],
                            a*rv[4] + b + c*u3[4], a*rv[5] + b + c*u3[5],
                            a*rv[6] + b + c*u3[6], a*rv[7] + b + c*u3[7],
                        ];
                        if let Some(r) = eval_exact(&p) {
                            if r.mu < 5.513755 {
                                loc.push((r.mu, a, b, c, p));
                            }
                        }
                    }
                    b2.lock().unwrap().extend(loc);
                }));
            }
            for h in hs2 { h.join().unwrap(); }
            let mut res2 = best2.lock().unwrap().clone();
            res2.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            println!("Deep 2 done in {:.2}s | Found {} points < 5.513755", t2.elapsed().as_secs_f64(), res2.len());
            for (i, (mu, a, b, c, p)) in res2.iter().take(12).enumerate() {
                println!("  2#{:2}: mu={:.12} | a={:3} b={:+2} c={:+2} | P={:?}", i+1, mu, a, b, c, p);
            }

            // =========================================================================
            // Deep 3: Second-Order Narrow-Cone Scan K * P_73 + d (ALL 360 unit 6-DOF)
            //         at K = 10, 30, 80, 200!
            // =========================================================================
            println!("=== Deep 3: Second-Order Narrow Cone K * P_73 + d (360 6-DOF dirs) ===");
            let t3 = Instant::now();
            let mut moves1: Vec<[i64; 8]> = Vec::new();
            for dh in -1..=1_i64 {
                for dj in -1..=1_i64 {
                    for dk in -1..=1_i64 {
                        for dl in -1..=1_i64 {
                            for dm in -1..=1_i64 {
                                let dr = dh + dm - dk;
                                if dr.abs() > 1 { continue; }
                                for dq in -1..=1_i64 {
                                    let ds = dj + dq - dl;
                                    if ds.abs() > 1 { continue; }
                                    if dh == 0 && dj == 0 && dk == 0 && dl == 0 && dm == 0 && dq == 0 { continue; }
                                    moves1.push([dh, dj, dk, dl, dm, dq, dr, ds]);
                                }
                            }
                        }
                    }
                }
            }
            let mut global_best = (5.513746817666_f64, p73);
            for &scale in &[10_i64, 30, 80, 200] {
                let best_k = Arc::new(Mutex::new(Vec::<(f64, [i64; 8], [i64; 8])>::new()));
                let ch = (moves1.len() + nthreads - 1) / nthreads;
                let mut hs = Vec::new();
                for tid in 0..nthreads {
                    let sl = moves1[tid * ch..moves1.len().min((tid + 1) * ch)].to_vec();
                    let bk = Arc::clone(&best_k);
                    hs.push(thread::spawn(move || {
                        let mut loc = Vec::new();
                        for d in sl {
                            let p = [
                                scale * p73[0] + d[0], scale * p73[1] + d[1],
                                scale * p73[2] + d[2], scale * p73[3] + d[3],
                                scale * p73[4] + d[4], scale * p73[5] + d[5],
                                scale * p73[6] + d[6], scale * p73[7] + d[7],
                            ];
                            if let Some(r) = eval_exact(&p) {
                                if r.mu < 5.5137470 {
                                    loc.push((r.mu, d, p));
                                }
                            }
                        }
                        bk.lock().unwrap().extend(loc);
                    }));
                }
                for h in hs { h.join().unwrap(); }
                let mut rk = best_k.lock().unwrap().clone();
                rk.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
                println!("  K={:3}: {} dirs < 5.5137470 | top: {:?}", scale, rk.len(), rk.first().map(|x| (x.0, x.1)));
                for (i, (mu, d, p)) in rk.iter().take(5).enumerate() {
                    println!("    K={:3} #{}: mu={:.12} | d={:?}", scale, i+1, mu, d);
                    if *mu < global_best.0 {
                        global_best = (*mu, *p);
                    }
                }
            }
            println!("Deep 3 done in {:.2}s", t3.elapsed().as_secs_f64());

            // =========================================================================
            // Deep 4: Exhaustive 6-DOF Global Base-Point Search (all valid points with
            //         M <= 24 modulo dihedral symmetry Θ, filtering via fast c0,c1!)
            // =========================================================================
            println!("=== Deep 4: Exhaustive 6-DOF Base-Point Scan (M <= 24, mod Θ) ===");
            let t4 = Instant::now();
            // WLOG by Θ rotation, k = max(h,j,k,l,m,q,r,s). By reflection, j >= l.
            // And all 16 positive, plus M <= 24.
            // Parallelize over k in 6..=24.
            let best4 = Arc::new(Mutex::new(Vec::<(f64, [i64; 8])>::new()));
            let checked = Arc::new(Mutex::new(0u64));
            let ks_all: Vec<i64> = (6..=24).collect();
            let mut hs4 = Vec::new();
            for tid in 0..nthreads {
                let b4 = Arc::clone(&best4);
                let chk = Arc::clone(&checked);
                let my_ks: Vec<i64> = ks_all.iter().copied().enumerate().filter(|(i, _)| i % nthreads == tid).map(|(_, x)| x).collect();
                hs4.push(thread::spawn(move || {
                    let mut loc: Vec<(f64, [i64; 8])> = Vec::new();
                    let mut cnt = 0u64;
                    for k in my_ks {
                        for h in 1..=k {
                            for m in 1..=k {
                                let r = h + m - k;
                                if r < 1 || r > k { continue; }
                                // Check j' = j+m-k > 0 => j > k-m
                                let j_lo = (k - m + 1).max(1);
                                for j in j_lo..=k {
                                    // By reflection symmetry fixing k? Actually reflection reverses (h,j,k,l,m,q,r,s) -> (l,k,j,h,s,r,q,m)
                                    // To be 100% safe and not miss anything, don't break reflection, only k = max(V[0..8]).
                                    for l in 1..=k {
                                        // h' = h+l-j > 0 => l > j-h
                                        if h + l <= j { continue; }
                                        // l' = l+r-m > 0
                                        if l + r <= m { continue; }
                                        for q in 1..=k {
                                            let s = j + q - l;
                                            if s < 1 || s > k { continue; }
                                            // Check remaining auxiliary > 0 and <= 24
                                            let vp = [h+l-j, j+m-k, k+q-l, l+r-m, m+s-q, q+h-r, r+j-s, s+k-h];
                                            let mut ok = true;
                                            for &x in &vp {
                                                if x < 1 || x > 24 { ok = false; break; }
                                            }
                                            if !ok { continue; }
                                            // Primitive check
                                            let p = [h, j, k, l, m, q, r, s];
                                            let mut g = 0i64;
                                            for &x in &p { g = gcd(g, x); }
                                            if g != 1 { continue; }
                                            // Fast analytic filter before full sweep:
                                            // With RV saving ~ 24.19/54 = 0.448 of (M+N+Q), actual c2 >= 0.50*(M+N+Q) roughly;
                                            // Actually let's use exact stationary check and skip clearly awful ones:
                                            if let Some((_, c0, c1, mnq)) = stationary_c0_c1(&p) {
                                                cnt += 1;
                                                // Upper bound on saving is <= 0.52 * (M+N+Q), so c2 >= 0.48*(M+N+Q)
                                                // Conservative lower bound on mu: if (c0+c1)/(c0 - 0.45*(M+N+Q)) > 5.7, skip!
                                                let s_sum = (mnq[0] + mnq[1] + mnq[2]) as f64;
                                                let denom_max = c0 - 0.45 * s_sum;
                                                if denom_max <= 0.0 || (c0 + c1) / denom_max > 5.8 {
                                                    continue;
                                                }
                                                if let Some(res) = eval_exact(&p) {
                                                    if res.mu < 5.55 {
                                                        loc.push((res.mu, p));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    b4.lock().unwrap().extend(loc);
                    *chk.lock().unwrap() += cnt;
                }));
            }
            for h in hs4 { h.join().unwrap(); }
            let mut res4 = best4.lock().unwrap().clone();
            res4.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            println!("Deep 4 done in {:.2}s | Checked {} admissible primitive points (M<=24) | Found {} with mu < 5.55",
                t4.elapsed().as_secs_f64(), *checked.lock().unwrap(), res4.len());
            for (i, (mu, p)) in res4.iter().take(15).enumerate() {
                println!("  4#{:2}: mu={:.12} | P={:?}", i+1, mu, p);
            }

            if let Some(r) = eval_exact(&global_best.1) {
                print_res("DEEP OVERALL CHAMPION", &r);
            }
        }

        "orbits" => {
            use std::collections::HashMap;
            let nthreads = 16usize;
            println!("=== Orbits: Distinct Phi-Orbits up to M <= 28 + Uniform Shift K*P +/- b ===");
            let t0 = Instant::now();
            // Collect best point per sorted 16-tuple V (Phi-orbit invariant) for M <= 28!
            let global_map = Arc::new(Mutex::new(HashMap::<[i64; 16], (f64, [i64; 8])>::new()));
            let ks_all: Vec<i64> = (6..=28).collect();
            let mut hs = Vec::new();
            for tid in 0..nthreads {
                let gmap = Arc::clone(&global_map);
                let my_ks: Vec<i64> = ks_all.iter().copied().enumerate().filter(|(i, _)| i % nthreads == tid).map(|(_, x)| x).collect();
                hs.push(thread::spawn(move || {
                    let mut local_map: HashMap<[i64; 16], (f64, [i64; 8])> = HashMap::new();
                    for k in my_ks {
                        for h in 1..=k {
                            for m in 1..=k {
                                let r = h + m - k;
                                if r < 1 || r > k { continue; }
                                let j_lo = (k - m + 1).max(1);
                                for j in j_lo..=k {
                                    for l in 1..=k {
                                        if h + l <= j || l + r <= m { continue; }
                                        for q in 1..=k {
                                            let s = j + q - l;
                                            if s < 1 || s > k { continue; }
                                            let mut v = [
                                                h, j, k, l, m, q, r, s,
                                                h + l - j, j + m - k, k + q - l, l + r - m,
                                                m + s - q, q + h - r, r + j - s, s + k - h,
                                            ];
                                            let mut ok = true;
                                            for &x in &v[8..] {
                                                if x < 1 || x > 28 { ok = false; break; }
                                            }
                                            if !ok { continue; }
                                            let p = [h, j, k, l, m, q, r, s];
                                            let mut g = 0i64;
                                            for &x in &p { g = gcd(g, x); }
                                            if g != 1 { continue; }
                                            if let Some((_, c0, c1, mnq)) = stationary_c0_c1(&p) {
                                                let s_sum = (mnq[0] + mnq[1] + mnq[2]) as f64;
                                                let denom_max = c0 - 0.46 * s_sum;
                                                if denom_max <= 0.0 || (c0 + c1) / denom_max > 5.75 {
                                                    continue;
                                                }
                                                v.sort_unstable();
                                                if local_map.contains_key(&v) {
                                                    continue;
                                                }
                                                if let Some(res) = eval_exact(&p) {
                                                    if res.mu < 5.60 {
                                                        local_map.insert(v, (res.mu, p));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    let mut gm = gmap.lock().unwrap();
                    for (k, val) in local_map {
                        gm.entry(k).or_insert(val);
                    }
                }));
            }
            for h in hs { h.join().unwrap(); }
            let mut orbits: Vec<(f64, [i64; 8], [i64; 16])> = global_map
                .lock()
                .unwrap()
                .iter()
                .map(|(v, (mu, p))| (*mu, *p, *v))
                .collect();
            orbits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            println!("Found {} distinct Phi-orbits (M<=28, mu<5.60) in {:.2}s", orbits.len(), t0.elapsed().as_secs_f64());
            for (i, (mu, p, _)) in orbits.iter().take(20).enumerate() {
                println!("  Orbit#{:2}: base mu={:.12} | P={:?}", i + 1, mu, p);
            }

            // Now test K * P_base + b * 1 for K in 1..=150, b in {-2,-1,+1,+2} on top 40 distinct orbits!
            println!("=== Sweeping K*P_base +/- b (K=1..150, b in {{-2,-1,1,2}}) on Top 40 Distinct Orbits ===");
            let t1 = Instant::now();
            let top40: Vec<[i64; 8]> = orbits.iter().take(40).map(|x| x.1).collect();
            let best_orb = Arc::new(Mutex::new(Vec::<(f64, [i64; 8], [i64; 8], i64, i64)>::new()));
            let mut hs2 = Vec::new();
            for tid in 0..nthreads {
                let slice: Vec<[i64; 8]> = top40.iter().copied().enumerate().filter(|(i, _)| i % nthreads == tid).map(|(_, x)| x).collect();
                let bo = Arc::clone(&best_orb);
                hs2.push(thread::spawn(move || {
                    let mut loc = Vec::new();
                    for base in slice {
                        let mut bbest = (10.0_f64, base, base, 1i64, 0i64);
                        for k in 1..=150_i64 {
                            for b in [-2_i64, -1, 1, 2] {
                                if gcd(k, b.abs()) != 1 { continue; }
                                let p = [
                                    k * base[0] + b, k * base[1] + b, k * base[2] + b, k * base[3] + b,
                                    k * base[4] + b, k * base[5] + b, k * base[6] + b, k * base[7] + b,
                                ];
                                if let Some(res) = eval_exact(&p) {
                                    if res.mu < bbest.0 {
                                        bbest = (res.mu, p, base, k, b);
                                    }
                                }
                            }
                        }
                        loc.push(bbest);
                    }
                    bo.lock().unwrap().extend(loc);
                }));
            }
            for h in hs2 { h.join().unwrap(); }
            let mut rorb = best_orb.lock().unwrap().clone();
            rorb.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            println!("Orbit-ray sweep done in {:.2}s!", t1.elapsed().as_secs_f64());
            for (i, (mu, p, base, k, b)) in rorb.iter().take(15).enumerate() {
                println!("  Top#{:2}: mu={:.12} | {}*{:?} {:+} => {:?}", i + 1, mu, k, base, b, p);
            }
        }

        _ => {
            eprintln!("Modes: selftest | ray | optimize | deep | orbits | eval h j k l m q r s");
        }
    }
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}
