//! expr_hunter — 16-core, low-memory (<40 MB RAM) Meet-in-the-Middle (MITM)
//! symbolic closed-form hunter searching >10^12 depth-4 symbolic expressions
//! over 40+ fundamental mathematical constants (pi, e, phi, sqrt(2), sqrt(3),
//! cbrt(3), ln(2), ln(phi), gamma, Catalan G, Gieseking V_tet, Glaisher A,
//! Gamma(1/4), Gamma(1/3), Gamma(1/6), eta(i), e^pi, e^(pi*sqrt(3)), AGM(1,sqrt2))
//! targeting zeta(3) and 14 nonlinear transforms of zeta(3).

use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

const ZETA3: f64 = 1.2020569031595942854;
const PI: f64 = std::f64::consts::PI;
const E: f64 = std::f64::consts::E;
const GAMMA_EULER: f64 = 0.57721566490153286061;
const CATALAN_G: f64 = 0.91596559417721901505;
const GIESEKING_V: f64 = 1.01494160640965362502;
const GLAISHER_A: f64 = 1.28242712910062263688;
const GAMMA_1_4: f64 = 3.62560990822190831193;
const GAMMA_1_3: f64 = 2.67893853470774763366;
const GAMMA_1_6: f64 = 5.56631600178023520425;
const ETA_I: f64 = 0.76822542232605663413;
const ETA_IS2: f64 = 0.69014833519717087855;
const ETA_IS3: f64 = 0.63554224055090801681;
const AGM_1_SQRT2: f64 = 1.19814023473559220744;
const OMEGA_W1: f64 = 0.56714329040978387300;
const PLASTIC_RHO: f64 = 1.32471795724474602596;
const KHINCHIN_K0: f64 = 2.68545200106530644531;

#[derive(Clone)]
struct D1Node {
    val: f64,
    cost: u8,
    expr: String,
}

#[derive(Clone, Copy)]
struct D2Node {
    val: f64,
    left: u16,
    right: u16,
    op: u8, // 0=D1-only, 1=+, 2=-, 3=*, 4=/, 5=^
    cost: u8,
}

fn format_d2(node: &D2Node, d1: &[D1Node]) -> String {
    if node.op == 0 {
        return d1[node.left as usize].expr.clone();
    }
    let l = &d1[node.left as usize].expr;
    let r = &d1[node.right as usize].expr;
    let sym = match node.op {
        1 => "+",
        2 => "-",
        3 => "*",
        4 => "/",
        5 => "**",
        _ => "?",
    };
    format!("({} {} {})", l, sym, r)
}

fn build_d1() -> Vec<D1Node> {
    let sqrt2 = 2.0_f64.sqrt();
    let sqrt3 = 3.0_f64.sqrt();
    let sqrt5 = 5.0_f64.sqrt();
    let phi = 0.5 * (1.0 + sqrt5);
    let cbrt2 = 2.0_f64.cbrt();
    let cbrt3 = 3.0_f64.cbrt();
    let cbrt5 = 5.0_f64.cbrt();

    let atoms: Vec<(f64, u8, &str)> = vec![
        (1.0, 1, "1"),
        (2.0, 1, "2"),
        (3.0, 1, "3"),
        (4.0, 2, "4"),
        (5.0, 2, "5"),
        (6.0, 2, "6"),
        (7.0, 2, "7"),
        (8.0, 2, "8"),
        (9.0, 2, "9"),
        (10.0, 2, "10"),
        (12.0, 2, "12"),
        (15.0, 3, "15"),
        (16.0, 2, "16"),
        (18.0, 3, "18"),
        (24.0, 3, "24"),
        (30.0, 3, "30"),
        (32.0, 3, "32"),
        (36.0, 3, "36"),
        (45.0, 3, "45"),
        (60.0, 3, "60"),
        (90.0, 3, "90"),
        (180.0, 3, "180"),
        (PI, 1, "pi"),
        (E, 2, "e"),
        (phi, 2, "phi"),
        (sqrt2, 2, "sqrt(2)"),
        (sqrt3, 2, "sqrt(3)"),
        (sqrt5, 2, "sqrt(5)"),
        (cbrt2, 3, "2**(1/3)"),
        (cbrt3, 3, "3**(1/3)"),
        (cbrt5, 3, "5**(1/3)"),
        (2.0_f64.ln(), 2, "ln(2)"),
        (3.0_f64.ln(), 2, "ln(3)"),
        (5.0_f64.ln(), 3, "ln(5)"),
        (phi.ln(), 2, "ln(phi)"),
        (PI.ln(), 3, "ln(pi)"),
        (GAMMA_EULER, 2, "gamma"),
        (CATALAN_G, 2, "G"),
        (GIESEKING_V, 3, "V_tet"),
        (GLAISHER_A, 3, "A"),
        (GAMMA_1_4, 3, "Gamma14"),
        (GAMMA_1_3, 3, "Gamma13"),
        (GAMMA_1_6, 3, "Gamma16"),
        (ETA_I, 3, "eta_i"),
        (ETA_IS2, 3, "eta_is2"),
        (ETA_IS3, 3, "eta_is3"),
        (PI.exp(), 3, "exp(pi)"),
        ((PI * sqrt2).exp(), 4, "exp(pi*sqrt(2))"),
        ((PI * sqrt3).exp(), 4, "exp(pi*sqrt(3))"),
        (AGM_1_SQRT2, 3, "AGM1s2"),
        (1.0 + sqrt2, 2, "silver"),
        ((1.0 + sqrt2).ln(), 2, "ln(1+sqrt(2))"),
        (6.0_f64.sqrt(), 3, "sqrt(6)"),
        (7.0_f64.sqrt(), 3, "sqrt(7)"),
        (OMEGA_W1, 3, "Omega"),
        (PLASTIC_RHO, 3, "plastic"),
        (KHINCHIN_K0, 3, "Khinchin"),
    ];

    let mut raw: Vec<D1Node> = Vec::new();
    for &(v, c, name) in &atoms {
        raw.push(D1Node { val: v, cost: c, expr: name.to_string() });
        // Unary ops on core atoms (cost <= 3)
        if c <= 3 {
            let unaries: Vec<(f64, u8, String)> = vec![
                (v * v, c + 1, format!("({})**2", name)),
                (v * v * v, c + 2, format!("({})**3", name)),
                (v.sqrt(), c + 2, format!("sqrt({})", name)),
                (v.cbrt(), c + 3, format!("cbrt({})", name)),
                (1.0 / v, c + 2, format!("(1/{})", name)),
                (v.ln(), c + 2, format!("ln({})", name)),
                (if v <= 15.0 { v.exp() } else { f64::NAN }, c + 2, format!("exp({})", name)),
                (v.sin(), c + 3, format!("sin({})", name)),
                (v.cos(), c + 3, format!("cos({})", name)),
                (v.tan(), c + 3, format!("tan({})", name)),
                (if v <= 12.0 { v.sinh() } else { f64::NAN }, c + 3, format!("sinh({})", name)),
                (if v <= 12.0 { v.cosh() } else { f64::NAN }, c + 3, format!("cosh({})", name)),
                (v.tanh(), c + 3, format!("tanh({})", name)),
                (v.atan(), c + 3, format!("atan({})", name)),
                (v.asinh(), c + 3, format!("asinh({})", name)),
            ];
            for (uv, uc, uexpr) in unaries {
                if uv.is_finite() && uv.abs() > 1e-10 && uv.abs() < 1e9 {
                    raw.push(D1Node { val: uv, cost: uc, expr: uexpr });
                }
            }
        }
    }

    // Deduplicate D1 nodes within 1e-13 keeping lowest cost
    raw.sort_by(|a, b| a.val.partial_cmp(&b.val).unwrap());
    let mut dedup: Vec<D1Node> = Vec::new();
    for node in raw {
        if let Some(last) = dedup.last_mut() {
            let tol = 1e-13 * node.val.abs().max(1.0);
            if (node.val - last.val).abs() < tol {
                if node.cost < last.cost || (node.cost == last.cost && node.expr.len() < last.expr.len()) {
                    *last = node;
                }
                continue;
            }
        }
        dedup.push(node);
    }
    dedup
}

fn build_d2(d1: &[D1Node]) -> Vec<D2Node> {
    // Expand D1 core (cost <= 6) up to pair cost <= 10 (~1.2M unique nodes -> ~20 MB RAM!)
    let core_idx: Vec<usize> = d1
        .iter()
        .enumerate()
        .filter(|(_, n)| n.cost <= 6)
        .map(|(i, _)| i)
        .collect();
    let pow_exp_idx: Vec<usize> = d1
        .iter()
        .enumerate()
        .filter(|(_, n)| n.cost <= 4 && n.val.abs() <= 5.0 && n.val != 0.0 && n.val != 1.0)
        .map(|(i, _)| i)
        .collect();

    println!(
        "D1 total: {} | D1 core (cost<=6): {} | pow exponents: {}",
        d1.len(),
        core_idx.len(),
        pow_exp_idx.len()
    );

    let mut raw: Vec<D2Node> = Vec::with_capacity(1_800_000);
    // Include all D1 nodes directly
    for (i, n) in d1.iter().enumerate() {
        raw.push(D2Node {
            val: n.val,
            left: i as u16,
            right: 0,
            op: 0,
            cost: n.cost,
        });
    }

    let push_val = |raw: &mut Vec<D2Node>, v: f64, l: usize, r: usize, op: u8, c: u8| {
        if v.is_finite() && v.abs() > 1e-9 && v.abs() < 1e7 {
            raw.push(D2Node {
                val: v,
                left: l as u16,
                right: r as u16,
                op,
                cost: c,
            });
        }
    };

    for (a_pos, &i) in core_idx.iter().enumerate() {
        let vi = d1[i].val;
        let ci = d1[i].cost;
        for &j in &core_idx[a_pos..] {
            let vj = d1[j].val;
            let cj = d1[j].cost;
            let c = ci.saturating_add(cj).saturating_add(1);
            if c > 10 {
                continue;
            }
            // Commutative: +, *
            push_val(&mut raw, vi + vj, i, j, 1, c);
            push_val(&mut raw, vi * vj, i, j, 3, c);
            // Non-commutative: -, / (both orders)
            push_val(&mut raw, vi - vj, i, j, 2, c);
            if i != j {
                push_val(&mut raw, vj - vi, j, i, 2, c);
            }
            if vj.abs() > 1e-8 {
                push_val(&mut raw, vi / vj, i, j, 4, c);
            }
            if i != j && vi.abs() > 1e-8 {
                push_val(&mut raw, vj / vi, j, i, 4, c);
            }
        }
        // Powers vi ** vj for positive vi and clean exponents vj
        if vi > 1e-6 && vi < 1e4 && (vi - 1.0).abs() > 1e-9 && ci <= 4 {
            for &j in &pow_exp_idx {
                let vj = d1[j].val;
                let c = ci.saturating_add(d1[j].cost).saturating_add(2);
                if c <= 10 {
                    let pv = vi.powf(vj);
                    push_val(&mut raw, pv, i, j, 5, c);
                }
            }
        }
    }

    println!("Raw D2 candidates before dedup: {}", raw.len());
    raw.sort_unstable_by(|a, b| a.val.total_cmp(&b.val));
    let mut dedup: Vec<D2Node> = Vec::with_capacity(raw.len() / 2);
    for node in raw {
        if let Some(last) = dedup.last_mut() {
            let tol = 1e-13 * node.val.abs().max(1.0);
            if (node.val - last.val).abs() < tol {
                if node.cost < last.cost {
                    *last = node;
                }
                continue;
            }
        }
        dedup.push(node);
    }
    dedup.shrink_to_fit();
    dedup
}

#[derive(Clone, Debug)]
struct Hit {
    target_name: &'static str,
    expr_str: String,
    rel_err: f64,
    total_cost: u8,
}

#[inline(always)]
fn lookup_match(table: &[D2Node], want: f64, rel_tol: f64) -> Option<usize> {
    if !want.is_finite() || want.abs() < 1e-9 || want.abs() > 1e7 {
        return None;
    }
    let idx = table.partition_point(|n| n.val < want);
    let abs_tol = rel_tol * want.abs().max(1e-6);
    let mut best_i: Option<usize> = None;
    let mut best_d = abs_tol;
    if idx < table.len() {
        let d = (table[idx].val - want).abs();
        if d < best_d {
            best_d = d;
            best_i = Some(idx);
        }
    }
    if idx > 0 {
        let d = (table[idx - 1].val - want).abs();
        if d < best_d {
            best_i = Some(idx - 1);
        }
    }
    best_i
}

fn main() {
    let t0 = Instant::now();
    let d1 = build_d1();
    let d2 = Arc::new(build_d2(&d1));
    let d1 = Arc::new(d1);
    let n = d2.len();
    println!(
        "Deduplicated D2 table: {} nodes (~{:.2} MB RAM) built in {:.2}s",
        n,
        (n * std::mem::size_of::<D2Node>()) as f64 / 1e6,
        t0.elapsed().as_secs_f64()
    );

    // Self-test: verify MITM instantly finds zeta(2) = pi^2/6 and beta(3) = pi^3/32
    for &(sname, sval) in &[("SELF-TEST pi^2/6", PI * PI / 6.0), ("SELF-TEST pi^3/32", PI * PI * PI / 32.0)] {
        if let Some(idx) = lookup_match(&d2, sval, 1e-13) {
            println!("[PASS] {} matched directly in D2: {} (err={:.2e})", sname, format_d2(&d2[idx], &d1), (d2[idx].val - sval).abs());
        } else {
            println!("[WARN] {} not direct in D2", sname);
        }
    }

    // Self-test #2: verify MITM finds genuine depth-4 identity Li_2(phi^-2) = pi^2/15 - ln(phi)^2!
    let phi = 0.5 * (1.0 + 5.0_f64.sqrt());
    let li2_phi_m2 = PI * PI / 15.0 - phi.ln() * phi.ln();
    let mut st_found = false;
    for l in d2.iter() {
        if let Some(j) = lookup_match(&d2, l.val - li2_phi_m2, 1e-14) {
            println!(
                "[PASS] SELF-TEST depth-4 Li_2(phi^-2) found via MITM: {} - {} (err={:.2e})",
                format_d2(l, &d1),
                format_d2(&d2[j], &d1),
                ((l.val - d2[j].val) - li2_phi_m2).abs()
            );
            st_found = true;
            break;
        }
    }
    assert!(st_found, "Depth-4 self-test failed!");

    let targets: Vec<(&'static str, f64)> = vec![
        ("zeta(3)", ZETA3),
        ("zeta(3)/pi^3", ZETA3 / (PI * PI * PI)),
        ("zeta(3)/pi^2", ZETA3 / (PI * PI)),
        ("zeta(3)/pi", ZETA3 / PI),
        ("1/zeta(3)", 1.0 / ZETA3),
        ("ln(zeta(3))", ZETA3.ln()),
        ("exp(zeta(3))", ZETA3.exp()),
        ("zeta(3)^2", ZETA3 * ZETA3),
        ("zeta(3)^3", ZETA3 * ZETA3 * ZETA3),
        ("sqrt(zeta(3))", ZETA3.sqrt()),
        ("sin(zeta(3))", ZETA3.sin()),
        ("cos(zeta(3))", ZETA3.cos()),
        ("tan(zeta(3))", ZETA3.tan()),
        ("sinh(zeta(3))", ZETA3.sinh()),
        ("cosh(zeta(3))", ZETA3.cosh()),
    ];

    let approx_space = (n as f64) * (n as f64) * 5.0 * (targets.len() as f64);
    println!(
        "=== Launching 16-Core MITM Search over ~{:.2e} Depth-4 Symbolic Expressions ===",
        approx_space
    );

    let t_search = Instant::now();
    let nthreads = 16usize;
    let hits = Arc::new(Mutex::new(Vec::<Hit>::new()));
    let chunk = (n + nthreads - 1) / nthreads;
    let mut handles = Vec::new();

    // Sweep ALL 3.14e13 D2 x D2 pairs at tol 2e-12 to surface both any exact hit AND the tightest symbolic near-misses!
    let rel_tol = 2.0e-12_f64;
    let max_total_cost = 22u8;

    for tid in 0..nthreads {
        let start = tid * chunk;
        let end = n.min((tid + 1) * chunk);
        let d2_ref = Arc::clone(&d2);
        let d1_ref = Arc::clone(&d1);
        let hits_ref = Arc::clone(&hits);
        let tgts = targets.clone();
        handles.push(thread::spawn(move || {
            let mut local_hits: Vec<Hit> = Vec::new();
            for i in start..end {
                let l = d2_ref[i];
                let lv = l.val;
                if l.cost > max_total_cost - 2 {
                    continue;
                }
                let ln_lv = if lv > 0.0 && (lv - 1.0).abs() > 1e-6 { Some(lv.ln()) } else { None };
                for &(tname, tv) in &tgts {
                    // 1. L + R = T
                    if let Some(j) = lookup_match(&d2_ref, tv - lv, rel_tol) {
                        if j >= i {
                            let r = d2_ref[j];
                            let tc = l.cost.saturating_add(r.cost);
                            if tc <= max_total_cost {
                                let v = lv + r.val;
                                let re = (v - tv).abs() / tv.abs().max(1.0);
                                if re < rel_tol && local_hits.len() < 5000 {
                                    local_hits.push(Hit {
                                        target_name: tname,
                                        expr_str: format!("{} + {}", format_d2(&l, &d1_ref), format_d2(&r, &d1_ref)),
                                        rel_err: re,
                                        total_cost: tc,
                                    });
                                }
                            }
                        }
                    }
                    // 2. L - R = T
                    if let Some(j) = lookup_match(&d2_ref, lv - tv, rel_tol) {
                        let r = d2_ref[j];
                        let tc = l.cost.saturating_add(r.cost);
                        if tc <= max_total_cost {
                            let v = lv - r.val;
                            let re = (v - tv).abs() / tv.abs().max(1.0);
                            if re < rel_tol && local_hits.len() < 5000 {
                                local_hits.push(Hit {
                                    target_name: tname,
                                    expr_str: format!("{} - {}", format_d2(&l, &d1_ref), format_d2(&r, &d1_ref)),
                                    rel_err: re,
                                    total_cost: tc,
                                });
                            }
                        }
                    }
                    // 3. L * R = T
                    if let Some(j) = lookup_match(&d2_ref, tv / lv, rel_tol) {
                        if j >= i {
                            let r = d2_ref[j];
                            let tc = l.cost.saturating_add(r.cost);
                            if tc <= max_total_cost {
                                let v = lv * r.val;
                                let re = (v - tv).abs() / tv.abs().max(1.0);
                                if re < rel_tol && local_hits.len() < 5000 {
                                    local_hits.push(Hit {
                                        target_name: tname,
                                        expr_str: format!("{} * {}", format_d2(&l, &d1_ref), format_d2(&r, &d1_ref)),
                                        rel_err: re,
                                        total_cost: tc,
                                    });
                                }
                            }
                        }
                    }
                    // 4. L / R = T
                    if let Some(j) = lookup_match(&d2_ref, lv / tv, rel_tol) {
                        let r = d2_ref[j];
                        let tc = l.cost.saturating_add(r.cost);
                        if tc <= max_total_cost {
                            let v = lv / r.val;
                            let re = (v - tv).abs() / tv.abs().max(1.0);
                            if re < rel_tol && local_hits.len() < 5000 {
                                local_hits.push(Hit {
                                    target_name: tname,
                                    expr_str: format!("{} / {}", format_d2(&l, &d1_ref), format_d2(&r, &d1_ref)),
                                    rel_err: re,
                                    total_cost: tc,
                                });
                            }
                        }
                    }
                    // 5. L ^ R = T
                    if tv > 0.0 {
                        if let Some(ln_l) = ln_lv {
                            let want_r = tv.ln() / ln_l;
                            if want_r.abs() <= 6.0 {
                                if let Some(j) = lookup_match(&d2_ref, want_r, rel_tol) {
                                    let r = d2_ref[j];
                                    let tc = l.cost.saturating_add(r.cost).saturating_add(1);
                                    if tc <= max_total_cost {
                                        let v = lv.powf(r.val);
                                        let re = (v - tv).abs() / tv.abs().max(1.0);
                                        if re < rel_tol && local_hits.len() < 5000 {
                                            local_hits.push(Hit {
                                                target_name: tname,
                                                expr_str: format!("({}) ** ({})", format_d2(&l, &d1_ref), format_d2(&r, &d1_ref)),
                                                rel_err: re,
                                                total_cost: tc,
                                            });
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            hits_ref.lock().unwrap().extend(local_hits);
        }));
    }
    for h in handles {
        h.join().unwrap();
    }

    let mut all_hits = hits.lock().unwrap().clone();
    // Sort by (total_cost, rel_err) so simplest & tightest appear first!
    all_hits.sort_by(|a, b| {
        a.total_cost
            .cmp(&b.total_cost)
            .then_with(|| a.rel_err.partial_cmp(&b.rel_err).unwrap())
    });

    println!(
        "MITM search done in {:.2}s! Total candidates (rel_err < {:.1e}, cost <= {}): {}",
        t_search.elapsed().as_secs_f64(),
        rel_tol,
        max_total_cost,
        all_hits.len()
    );
    println!("\n=== TOP 40 SIMPLEST CANDIDATE HITS (MUST BE VERIFIED TO 120 DIGITS IN PYTHON!) ===");
    for (k, h) in all_hits.iter().take(40).enumerate() {
        println!(
            "HIT#{:2} [cost={:2}] {} == {} | err={:.2e}",
            k + 1,
            h.total_cost,
            h.target_name,
            h.expr_str,
            h.rel_err
        );
    }
    // Also print top 20 by pure smallest rel_err
    all_hits.sort_by(|a, b| a.rel_err.partial_cmp(&b.rel_err).unwrap());
    println!("\n=== TOP 20 TIGHTEST RESIDUAL CANDIDATES ===");
    for (k, h) in all_hits.iter().take(20).enumerate() {
        println!(
            "TIGHT#{:2} [cost={:2}] {} == {} | err={:.2e}",
            k + 1,
            h.total_cost,
            h.target_name,
            h.expr_str,
            h.rel_err
        );
    }
}
