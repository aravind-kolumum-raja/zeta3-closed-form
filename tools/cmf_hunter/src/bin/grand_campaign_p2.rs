// grand_campaign_p2.rs
// ====================
// Phase 2 of the 25-Experiment Grand Campaign (Exps 6, 7, 8, 9, 10, 21):
// Searches ALL previously untouched 2-block double-pole hypergeometric families:
//   - Exp 6:  Coincident zero starts (z1 = z2 = 0 => t0 = 0!) and coincident ends
//   - Exp 7:  Outer pole s1 = 0 starting at t = 0
//   - Exp 8:  Unbalanced / super-convergent degree excess Delta_def = (D1+D2) - sum(L) - sum(K) > 0
//   - Exp 9:  High-L / asymmetric unequal-L 3-block (L up to 60) + multi-start evolutionary optimizer
//   - Exp 10: Deep unequal-L 4-, 5-, 6-block double-pole families
//   - Exp 21: TWO-SIDED ZERO BLOCKS (right zeros + LEFT zeros at t < -p1*n, like Beukers/Zudilin!)
//
// Generalized rational function:
//   F_n(t) = [ prod_{u=1}^{m_R} (t - r_u n)_{L_u n} * prod_{v=1}^{m_L} (t + q_v n + 1)_{K_v n} ]
//            / [ (t + s1*n)_{D1*n + 1} * (t + s2*n)_{D2*n + 1} ]
//
// Uses 16-core parallel screening + exact Farey certification (3 primes/cell) +
// multi-start simulated annealing / coordinate descent directly on exact mu!

use std::collections::HashSet;
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct GenSpec {
    // Right zero blocks: (z_u, r_u) with 0 <= z_u < r_u, L_u = r_u - z_u
    r_blocks: Vec<(i32, i32)>,
    // Left zero blocks: (q_v, K_v) with q_v >= p1, K_v >= 1 -> zeros [-q_v - K_v, -q_v - 1]
    l_blocks: Vec<(i32, i32)>,
    // Nested poles: 0 <= s1 <= s2 <= p2 <= p1, D1 = p1 - s1, D2 = p2 - s2
    s1: i32,
    p1: i32,
    s2: i32,
    p2: i32,
}

impl GenSpec {
    fn d1(&self) -> i32 { self.p1 - self.s1 }
    fn d2(&self) -> i32 { self.p2 - self.s2 }
    fn e1(&self) -> i32 { self.p1 - self.s2 }
    fn e2(&self) -> i32 { self.p2 - self.s1 }
    fn deg_def(&self) -> i32 {
        let sum_r: i32 = self.r_blocks.iter().map(|&(z, r)| r - z).sum();
        let sum_l: i32 = self.l_blocks.iter().map(|&(_, k)| k).sum();
        (self.d1() + self.d2()) - (sum_r + sum_l)
    }
    fn canonical(&self) -> GenSpec {
        let mut s = self.clone();
        s.r_blocks.sort();
        s.l_blocks.sort();
        s
    }
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
    let rad = (1.0 + c[..d].iter().map(|x| x.abs()).fold(0.0_f64, f64::max)).min(120.0);
    let mut roots: Vec<Cpx> = (0..d).map(|k| {
        let ang = 2.0 * PI * (k as f64 + 0.13719) / (d as f64);
        Cpx::new(rad * ang.cos(), rad * ang.sin())
    }).collect();
    for _ in 0..350 {
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

#[derive(Clone, Debug)]
struct Saddle {
    sigma: f64,
    g: f64,
    t0: i32,
    m_max: i32,
    z_re: f64,
    z_im: f64,
}

fn eval_saddle(sp: &GenSpec) -> Option<Saddle> {
    if sp.r_blocks.len() < 2 { return None; }
    if sp.s1 < 0 || sp.s2 < sp.s1 || sp.p2 < sp.s2 || sp.p1 < sp.p2 { return None; }
    if sp.d1() <= 0 || sp.d2() <= 0 || sp.e1() <= 0 || sp.e2() <= 0 { return None; }
    if sp.deg_def() < 0 { return None; }
    for &(z, r) in &sp.r_blocks { if z < 0 || r <= z { return None; } }
    for &(q, k) in &sp.l_blocks { if q < sp.p1 || k <= 0 { return None; } }

    let mut zs: Vec<i32> = sp.r_blocks.iter().map(|&(z, _)| z).collect();
    zs.sort();
    let mut rs: Vec<i32> = sp.r_blocks.iter().map(|&(_, r)| r).collect();
    rs.sort();
    let t0 = zs[1];
    let t1 = rs[rs.len() - 2];
    if t1 <= t0 { return None; }

    // Saddle eq: prod_u(z-z_u) * prod_v(z+q_v+K_v) * (z+s1)(z+s2)
    //          = prod_u(z-r_u) * prod_v(z+q_v)     * (z+p1)(z+p2)
    let mut p_l = vec![1.0_f64];
    for &(z, _) in &sp.r_blocks { p_l = poly_mul(&p_l, &[-(z as f64), 1.0]); }
    for &(q, k) in &sp.l_blocks { p_l = poly_mul(&p_l, &[(q + k) as f64, 1.0]); }
    p_l = poly_mul(&p_l, &[sp.s1 as f64, 1.0]);
    p_l = poly_mul(&p_l, &[sp.s2 as f64, 1.0]);

    let mut p_r = vec![1.0_f64];
    for &(_, r) in &sp.r_blocks { p_r = poly_mul(&p_r, &[-(r as f64), 1.0]); }
    for &(q, _) in &sp.l_blocks { p_r = poly_mul(&p_r, &[q as f64, 1.0]); }
    p_r = poly_mul(&p_r, &[sp.p1 as f64, 1.0]);
    p_r = poly_mul(&p_r, &[sp.p2 as f64, 1.0]);

    let diff: Vec<f64> = p_l.iter().zip(p_r.iter()).map(|(a, b)| a - b).collect();
    let roots = durand_kerner(diff);

    let e1 = sp.e1() as f64;
    let e2 = sp.e2() as f64;
    let mut c_const = e1 * e1.ln() + e2 * e2.ln();
    for &(z, r) in &sp.r_blocks { let l = (r - z) as f64; c_const -= l * l.ln(); }
    for &(_, k) in &sp.l_blocks { let kf = k as f64; c_const -= kf * kf.ln(); }

    let mut sigma_best = f64::NEG_INFINITY;
    let mut err_best = f64::NEG_INFINITY;
    let mut z_re_best = 0.0;
    let mut z_im_best = 0.0;

    // Evaluate pole saddle(s) in (-p2, -s2)
    for z in &roots {
        if z.im.abs() < 1e-6 && z.re > -(sp.p2 as f64) + 1e-7 && z.re < -(sp.s2 as f64) - 1e-7 {
            let j = -z.re;
            let mut sig = c_const;
            let mut ok = true;
            for &(zu, ru) in &sp.r_blocks {
                let a = j + ru as f64; let b = j + zu as f64;
                if a <= 0.0 || b < 0.0 { ok = false; break; }
                sig += a * a.ln() - if b > 1e-14 { b * b.ln() } else { 0.0 };
            }
            for &(qv, kv) in &sp.l_blocks {
                let a = (qv + kv) as f64 - j; let b = qv as f64 - j;
                if a <= 0.0 || b < 0.0 { ok = false; break; }
                sig += a * a.ln() - if b > 1e-14 { b * b.ln() } else { 0.0 };
            }
            let d1a = j - sp.s1 as f64; let d1b = sp.p1 as f64 - j;
            let d2a = j - sp.s2 as f64; let d2b = sp.p2 as f64 - j;
            if !ok || d1a <= 0.0 || d1b <= 0.0 || d2a <= 0.0 || d2b <= 0.0 { continue; }
            sig -= d1a * d1a.ln() + d1b * d1b.ln() + d2a * d2a.ln() + d2b * d2b.ln();
            if sig > sigma_best { sigma_best = sig; }
        }
    }
    if !sigma_best.is_finite() { return None; }

    // Evaluate upper half-plane complex saddle(s) with Re(z) in (-s1, t1) — or any Im > 0 saddle!
    // Also check real positive saddle on (t1, inf) if deg_def > 0!
    // WAIT: When deg_def > 0, check if there is a REAL saddle on (max_r, inf) OR a complex saddle!
    let max_r = *rs.last().unwrap() as f64;
    for z in &roots {
        if z.im > 1e-6 {
            let mut val = Cpx::new(c_const, 0.0);
            for &(zu, ru) in &sp.r_blocks {
                val = val.add(Cpx::new(z.re - zu as f64, z.im).xlogx())
                         .sub(Cpx::new(z.re - ru as f64, z.im).xlogx());
            }
            for &(qv, kv) in &sp.l_blocks {
                val = val.add(Cpx::new(z.re + (qv + kv) as f64, z.im).xlogx())
                         .sub(Cpx::new(z.re + qv as f64, z.im).xlogx());
            }
            val = val.sub(Cpx::new(z.re + sp.p1 as f64, z.im).xlogx())
                     .add(Cpx::new(z.re + sp.s1 as f64, z.im).xlogx())
                     .sub(Cpx::new(z.re + sp.p2 as f64, z.im).xlogx())
                     .add(Cpx::new(z.re + sp.s2 as f64, z.im).xlogx());
            let err = val.re - 2.0 * PI * z.im;
            if err > err_best {
                err_best = err;
                z_re_best = z.re;
                z_im_best = z.im;
            }
        } else if z.im.abs() < 1e-6 && z.re > max_r + 1e-5 {
            // Real positive saddle in the tail sum t > max_r! (Occurs when deg_def > 0!)
            let x = z.re;
            let mut val = c_const;
            for &(zu, ru) in &sp.r_blocks {
                val += (x - zu as f64) * (x - zu as f64).ln() - (x - ru as f64) * (x - ru as f64).ln();
            }
            for &(qv, kv) in &sp.l_blocks {
                val += (x + (qv + kv) as f64) * (x + (qv + kv) as f64).ln() - (x + qv as f64) * (x + qv as f64).ln();
            }
            val -= (x + sp.p1 as f64) * (x + sp.p1 as f64).ln() - (x + sp.s1 as f64) * (x + sp.s1 as f64).ln();
            val -= (x + sp.p2 as f64) * (x + sp.p2 as f64).ln() - (x + sp.s2 as f64) * (x + sp.s2 as f64).ln();
            if val > err_best {
                err_best = val;
                z_re_best = x;
                z_im_best = 0.0;
            }
        }
    }
    if !err_best.is_finite() { return None; }
    let g = sigma_best - err_best;
    if g <= 0.0 { return None; }
    Some(Saddle { sigma: sigma_best, g, t0, m_max: sp.p1 + t0, z_re: z_re_best, z_im: z_im_best })
}

fn sieve_primes(limit: usize) -> Vec<u64> {
    let mut is_p = vec![true; limit + 1];
    is_p[0] = false; is_p[1] = false;
    let mut i = 2;
    while i * i <= limit {
        if is_p[i] { let mut j = i * i; while j <= limit { is_p[j] = false; j += i; } }
        i += 1;
    }
    (2..=limit).filter(|&x| is_p[x]).map(|x| x as u64).collect()
}

fn mod_inv(a: i128, m: i128) -> i128 {
    let (mut t, mut nt, mut r, mut nr) = (0i128, 1i128, m, a.rem_euclid(m));
    while nr != 0 {
        let q = r / nr;
        let tmp = t - q * nt; t = nt; nt = tmp;
        let tmp = r - q * nr; r = nr; nr = tmp;
    }
    t.rem_euclid(m)
}

#[inline]
fn vp_of_u64(mut x: u64, p: u64) -> (i32, u64) {
    let mut v = 0;
    while x % p == 0 { v += 1; x /= p; }
    (v, x)
}

// Compute exact p-adic valuation min(v_p(2 d_{Mn}^3 C_n Q_n), v_p(2 d_{Mn}^3 C_n P_{3,n}))
// using Z/p^K Z arithmetic (K = 5 for safety!).
// Valid whenever p > 2 and every factorial argument < p^2 (or even general via exact table!).
fn vp_exact_gen(sp: &GenSpec, n: usize, p: u64, t0: i32) -> i32 {
    let max_r = sp.r_blocks.iter().map(|&(_, r)| r).max().unwrap() as usize;
    let max_l = sp.l_blocks.iter().map(|&(q, k)| q + k).max().unwrap_or(0) as usize;
    let p1n = sp.p1 as usize * n;
    let p2n = sp.p2 as usize * n;
    let s1n = sp.s1 as usize * n;
    let s2n = sp.s2 as usize * n;
    let t0n = t0 as usize * n;
    let mn = p1n + t0n;
    let lim = (p1n + max_r * n).max(max_l * n).max(mn) + 5;

    let k_prec: u32 = if p <= 800 { 5 } else { 4 };
    let pk: i128 = (p as i128).pow(k_prec);
    let mut vfac = vec![0i32; lim + 1];
    let mut ufac = vec![1i128; lim + 1];
    for i in 1..=lim {
        let (vi, ui) = vp_of_u64(i as u64, p);
        vfac[i] = vfac[i - 1] + vi;
        ufac[i] = (ufac[i - 1] * (ui as i128)) % pk;
    }

    // Harmonic sums H_N^{(1)}, H_N^{(2)}, H_N^{(3)} scaled by p^1, p^2, p^3 (valid when N < p^2!)
    if (lim as u64) >= p * p { return -1; }
    let p128 = p as i128;
    let p2_128 = p128 * p128;
    let p3_128 = p2_128 * p128;
    let mut p_h1 = vec![0i128; lim + 1];
    let mut p2_h2 = vec![0i128; mn + 1];
    let mut p3_h3 = vec![0i128; mn + 1];
    for k in 1..=lim {
        let (vk, uk) = vp_of_u64(k as u64, p);
        let inv1 = mod_inv(uk as i128, pk);
        let add1 = if vk == 0 { (p128 * inv1) % pk } else { inv1 };
        p_h1[k] = (p_h1[k - 1] + add1) % pk;
        if k <= mn {
            let inv2 = (inv1 * inv1) % pk;
            let inv3 = (inv2 * inv1) % pk;
            let add2 = if vk == 0 { (p2_128 * inv2) % pk } else { inv2 };
            let add3 = if vk == 0 { (p3_128 * inv3) % pk } else { inv3 };
            p2_h2[k] = (p2_h2[k - 1] + add2) % pk;
            p3_h3[k] = (p3_h3[k - 1] + add3) % pk;
        }
    }

    // Normalization C_n = (E1*n)!(E2*n)! / [ prod (L_u*n)! * prod (K_v*n)! ]
    let e1n = sp.e1() as usize * n;
    let e2n = sp.e2() as usize * n;
    let mut v_cn = vfac[e1n] + vfac[e2n];
    let mut u_cn_num = (ufac[e1n] * ufac[e2n]) % pk;
    let mut u_cn_den = 1i128;
    for &(zu, ru) in &sp.r_blocks {
        let ln = (ru - zu) as usize * n;
        v_cn -= vfac[ln];
        u_cn_den = (u_cn_den * ufac[ln]) % pk;
    }
    for &(_, kv) in &sp.l_blocks {
        let kn = kv as usize * n;
        v_cn -= vfac[kn];
        u_cn_den = (u_cn_den * ufac[kn]) % pk;
    }
    let u_cn = (u_cn_num * mod_inv(u_cn_den, pk)) % pk;

    // Parity sign factor for double poles: (-1)^{sum L_u n + (s1+s2)n} -> constant, omit or include!
    // Wait: we MUST keep the EXACT relative sign between double poles and simple-pole wings!
    // Double pole j in [s2n, p2n]:
    //   sign_dbl = (-1)^{sum L_u*n + s1n + s2n}
    // Left wing j in [s1n, s2n):
    //   sign_lw(j) = (-1)^{sum L_u*n + j - s1n}  => relative to sign_dbl: (-1)^{j - s2n}
    // Right wing j in (p2n, p1n]:
    //   sign_rw(j) = (-1)^{sum L_u*n + j - s1n + D2*n + 1} => relative to sign_dbl: (-1)^{j - p2n + 1}
    let inv2 = mod_inv(2, pk);
    let d3_pow = if (mn as u64) >= p { 3 } else { 0 }; // v_p(d_{Mn}^3) since Mn < p^2

    // Accumulate d_{Mn}^3 C_n Q_n and d_{Mn}^3 C_n P_{3,n} with floating minimum exponent!
    let mut q_terms: Vec<(i32, i128)> = Vec::with_capacity(p1n - s1n + 2);
    let mut p_terms: Vec<(i32, i128)> = Vec::with_capacity((p1n - s1n + 2) * 2);

    for j in s1n..=p1n {
        let in2 = j >= s2n && j <= p2n;
        // Numerator factorial product
        let mut vn = v_cn;
        let mut un_top = u_cn;
        let mut un_bot = 1i128;
        for &(zu, ru) in &sp.r_blocks {
            let a = j + ru as usize * n;
            let b = j + zu as usize * n;
            vn += vfac[a] - vfac[b];
            un_top = (un_top * ufac[a]) % pk;
            un_bot = (un_bot * ufac[b]) % pk;
        }
        for &(qv, kv) in &sp.l_blocks {
            let a = (qv + kv) as usize * n - j;
            let b = qv as usize * n - j;
            vn += vfac[a] - vfac[b];
            un_top = (un_top * ufac[a]) % pk;
            un_bot = (un_bot * ufac[b]) % pk;
        }
        // Denominator block 1 omitting j: (j-s1n)!(p1n-j)!
        let a1 = j - s1n; let b1 = p1n - j;
        vn -= vfac[a1] + vfac[b1];
        un_bot = (un_bot * ((ufac[a1] * ufac[b1]) % pk)) % pk;

        if in2 {
            // Double pole: block 2 omitting j is (j-s2n)!(p2n-j)!
            let a2 = j - s2n; let b2 = p2n - j;
            vn -= vfac[a2] + vfac[b2];
            un_bot = (un_bot * ((ufac[a2] * ufac[b2]) % pk)) % pk;
            let u_aj = (un_top * mod_inv(un_bot, pk)) % pk;

            // Q_n term: d_{Mn}^3 * C_n a_j has valuation vn + d3_pow
            q_terms.push((vn + d3_pow, u_aj));

            // P3_n part 1: d_{Mn}^3 * C_n a_j * H_{j+t0n}^{(3)}
            // Note: p^3 * H^{(3)} = p3_h3[j+t0n], so d_{Mn}^3 * H^{(3)} = p^{d3_pow - 3} * p3_h3[j+t0n]!
            let hj = j + t0n;
            let h3_val = p3_h3[hj];
            if h3_val != 0 {
                p_terms.push((vn + d3_pow - 3, (u_aj * h3_val) % pk));
            }

            // P3_n part 2: (1/2) * d_{Mn}^3 * C_n b_j * H_{j+t0n}^{(2)}
            // p * (b_j / a_j) = sum_u(pH[j+zu]-pH[j+ru]) + sum_v(pH[qv+kv-j]-pH[qv-j]) + pH[a1]-pH[b1] + pH[a2]-pH[b2]
            let mut pld: i128 = 0;
            for &(zu, ru) in &sp.r_blocks {
                pld += p_h1[j + zu as usize * n] - p_h1[j + ru as usize * n];
            }
            for &(qv, kv) in &sp.l_blocks {
                pld += p_h1[(qv + kv) as usize * n - j] - p_h1[qv as usize * n - j];
            }
            pld += p_h1[a1] - p_h1[b1] + p_h1[a2] - p_h1[b2];
            pld = pld.rem_euclid(pk);
            let h2_val = p2_h2[hj];
            // (1/2) * a_j * (pld / p) * (h2_val / p^2) * p^{d3_pow} => exp = vn + d3_pow - 3
            let u_part2 = (((((u_aj * pld) % pk) * h2_val) % pk) * inv2) % pk;
            if u_part2 != 0 {
                p_terms.push((vn + d3_pow - 3, u_part2));
            }
        } else if j < s2n {
            // Left simple-pole wing: block 2 = (p2n-j)! / (s2n-j-1)!
            let top2 = p2n - j; let bot2 = s2n - j - 1;
            vn -= vfac[top2] - vfac[bot2];
            un_top = (un_top * ufac[bot2]) % pk;
            un_bot = (un_bot * ufac[top2]) % pk;
            let mut u_bj = (un_top * mod_inv(un_bot, pk)) % pk;
            if (s2n - j) % 2 == 1 { u_bj = (pk - u_bj) % pk; }
            let hj = j + t0n;
            let h2_val = p2_h2[hj];
            // (1/2) * b_j * (h2_val / p^2) * p^{d3_pow} => exp = vn + d3_pow - 2
            let u_part = (((u_bj * h2_val) % pk) * inv2) % pk;
            if u_part != 0 {
                p_terms.push((vn + d3_pow - 2, u_part));
            }
        } else {
            // Right simple-pole wing (j > p2n): block 2 = (-1)^{D2*n+1} * (j-s2n)! / (j-p2n-1)!
            let top2 = j - s2n; let bot2 = j - p2n - 1;
            vn -= vfac[top2] - vfac[bot2];
            un_top = (un_top * ufac[bot2]) % pk;
            un_bot = (un_bot * ufac[top2]) % pk;
            let mut u_bj = (un_top * mod_inv(un_bot, pk)) % pk;
            if (j - p2n + 1) % 2 == 1 { u_bj = (pk - u_bj) % pk; }
            let hj = j + t0n;
            let h2_val = p2_h2[hj];
            let u_part = (((u_bj * h2_val) % pk) * inv2) % pk;
            if u_part != 0 {
                p_terms.push((vn + d3_pow - 2, u_part));
            }
        }
    }

    let eval_vec_vp = |terms: &[(i32, i128)]| -> i32 {
        if terms.is_empty() { return 99; }
        let vmin = terms.iter().map(|&(v, _)| v).min().unwrap();
        let mut acc: i128 = 0;
        for &(v, u) in terms {
            let dv = (v - vmin) as u32;
            if dv < k_prec {
                acc = (acc + u * p128.pow(dv)) % pk;
            }
        }
        let mut ext = 0i32;
        if acc == 0 { return vmin + k_prec as i32; }
        while acc % p128 == 0 && ext < k_prec as i32 {
            ext += 1;
            acc /= p128;
        }
        vmin + ext
    };

    let vq = eval_vec_vp(&q_terms);
    let vp = eval_vec_vp(&p_terms);
    vq.min(vp)
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
    for d in 2..=n { for k in 1..d { if gcd_i32(k, d) == 1 { pts.push((k, d)); } } }
    pts.sort_by(|a, b| (a.0 as i64 * b.1 as i64).cmp(&(b.0 as i64 * a.1 as i64)));
    pts
}

fn max_arg_1(sp: &GenSpec, sr: &Saddle) -> u64 {
    let max_r = sp.r_blocks.iter().map(|&(_, r)| r).max().unwrap();
    let max_l = sp.l_blocks.iter().map(|&(q, k)| q + k).max().unwrap_or(0);
    ((sp.p1 + max_r).max(max_l).max(sr.m_max)) as u64
}

fn fast_mu(sp: &GenSpec, sr: &Saddle, primes: &[u64]) -> Option<(f64, f64)> {
    let ma1 = max_arg_1(sp, sr);
    if 2 * ma1 + 15 > 2400 { return None; }
    let bins = 40;
    let dw = (sr.m_max as f64 - 1.0) / (bins as f64);
    let mut d1 = 0.0;
    for bi in 0..bins {
        let w = 1.0 + (bi as f64 + 0.5) * dw;
        let u = 1.0 / w;
        let min_p = (((ma1 as f64) * u).ceil() as u64 + 6).max(sr.m_max as u64 + 5);
        let mut best: Option<(usize, u64, f64)> = None;
        for &p in primes {
            if p < min_p { continue; }
            if p > min_p + 180 || p > 2450 { break; }
            let n = (u * p as f64).round() as usize;
            if n == 0 || n >= p as usize { continue; }
            let err = ((n as f64) / (p as f64) - u).abs();
            if best.map_or(true, |b| err < b.2) {
                best = Some((n, p, err));
                if err < 0.002 * u { break; }
            }
        }
        let (n, p, _) = best?;
        let v = vp_exact_gen(sp, n, p, sr.t0);
        if v < 0 { return None; }
        d1 += (v as f64) * dw;
    }
    let p_per = *primes.iter().find(|&&q| q >= 2 * ma1 + 11 && q <= 2450)?;
    let mut d2 = 0.0;
    for bi in 0..16 {
        let a = bi as f64 / 16.0; let b = (bi + 1) as f64 / 16.0;
        let k = (((a + b) * 0.5) * p_per as f64).round() as usize;
        let k = k.max(1).min(p_per as usize - 1);
        let v = vp_exact_gen(sp, p_per as usize + k, p_per, sr.t0);
        if v > 0 { d2 += (v as f64) * (digamma_p1(b) - digamma_p1(a)); }
    }
    let dt = d1 + d2;
    let h = sr.sigma + (3 * sr.m_max) as f64 - dt;
    let mu = if sr.g > h + 1e-6 { sr.g / (sr.g - h) } else { 999.0 };
    Some((mu, dt))
}

fn exact_farey_mu(sp: &GenSpec, sr: &Saddle, primes: &[u64], n_primes: usize) -> Option<(f64, f64, f64, f64, usize)> {
    let ma1 = max_arg_1(sp, sr);
    let farey_n = (ma1 as i32 + 4).min(65);
    let fpts = farey_seq(farey_n);
    let inv_m = 1.0 / (sr.m_max as f64);
    let mut d1 = 0.0;
    let mut dis = 0usize;
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        if b <= inv_m + 1e-12 { continue; }
        let a_cl = a.max(inv_m);
        if b - a_cl < 1e-12 { continue; }
        let min_p = (((ma1 as f64) * b).ceil() as u64 + 6).max(sr.m_max as u64 + 6);
        let mut vals = Vec::new();
        for &p in primes {
            if p < min_p || p > 2950 { continue; }
            let pf = p as f64; let mrg = (b - a_cl) * 0.08;
            let lo = ((a_cl + mrg) * pf).ceil() as i64;
            let hi = ((b - mrg) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let n = ((lo + hi) / 2) as usize;
                if p * p > ma1 * (n as u64) + 5 {
                    let v = vp_exact_gen(sp, n, p, sr.t0);
                    if v < 0 { return None; }
                    vals.push(v);
                    if vals.len() >= n_primes { break; }
                }
            }
        }
        if vals.is_empty() { return None; }
        let vmin = *vals.iter().min().unwrap();
        let vmax = *vals.iter().max().unwrap();
        if vmin != vmax { dis += 1; }
        if vmin > 0 { d1 += (vmin as f64) * (1.0 / a_cl - 1.0 / b); }
    }
    let min_p_per = 2 * ma1 + 11;
    let mut d2 = 0.0;
    for i in 0..fpts.len() - 1 {
        let a = fpts[i].0 as f64 / fpts[i].1 as f64;
        let b = fpts[i + 1].0 as f64 / fpts[i + 1].1 as f64;
        let mut vals = Vec::new();
        for &p in primes {
            if p < min_p_per || p > 2950 { continue; }
            let pf = p as f64; let mrg = (b - a) * 0.08;
            let lo = ((a + mrg) * pf).ceil() as i64;
            let hi = ((b - mrg) * pf).floor() as i64;
            if lo <= hi && lo >= 1 {
                let k = ((lo + hi) / 2) as usize;
                let n = p as usize + k;
                if p * p > ma1 * (n as u64) + 5 {
                    let v = vp_exact_gen(sp, n, p, sr.t0);
                    if v < 0 { return None; }
                    vals.push(v);
                    if v == 0 || vals.len() >= 2 { break; }
                }
            }
        }
        if vals.is_empty() { return None; }
        let vmin = *vals.iter().min().unwrap();
        if vmin > 0 { d2 += (vmin as f64) * (digamma_p1(b) - digamma_p1(a)); }
    }
    let dt = d1 + d2;
    let h = sr.sigma + (3 * sr.m_max) as f64 - dt;
    let mu = if sr.g > h + 1e-6 { sr.g / (sr.g - h) } else { 999.0 };
    Some((mu, dt, d1, d2, dis))
}

// Fast PRNG (xorshift64)
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        self.0 = x; x
    }
    fn range(&mut self, lo: i32, hi: i32) -> i32 {
        if hi <= lo { return lo; }
        lo + (self.next() % ((hi - lo + 1) as u64)) as i32
    }
    fn f01(&mut self) -> f64 { (self.next() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64)) }
}

fn mutate_spec(sp: &GenSpec, rng: &mut Rng, allow_unbal: bool, allow_left: bool) -> GenSpec {
    let mut s = sp.clone();
    let n_muts = if rng.f01() < 0.65 { 1 } else { 2 };
    for _ in 0..n_muts {
        let choice = rng.range(0, 7);
        match choice {
            0 | 1 => {
                // Perturb a right block endpoint z or r (keeping deg_def valid if balanced)
                let idx = rng.range(0, s.r_blocks.len() as i32 - 1) as usize;
                let step = if rng.f01() < 0.5 { 1 } else { -1 };
                if rng.f01() < 0.5 {
                    // Shift whole block (z, r) -> (z+step, r+step) (preserves L!)
                    if s.r_blocks[idx].0 + step >= 0 {
                        s.r_blocks[idx].0 += step;
                        s.r_blocks[idx].1 += step;
                    }
                } else {
                    // Change length L_idx by step, compensate in another block or p1/p2
                    if rng.f01() < 0.5 { s.r_blocks[idx].1 += step; } else { s.r_blocks[idx].0 -= step; }
                    if !allow_unbal {
                        // Compensate in another block or pole
                        let c2 = rng.range(0, 2);
                        if c2 == 0 { s.p1 += step; }
                        else if c2 == 1 { s.p2 += step; }
                        else {
                            let j = rng.range(0, s.r_blocks.len() as i32 - 1) as usize;
                            if j != idx { s.r_blocks[j].1 -= step; }
                            else { s.p1 += step; }
                        }
                    }
                }
            }
            2 | 3 => {
                // Perturb pole endpoints s1, s2, p1, p2
                let step = if rng.f01() < 0.5 { 1 } else { -1 };
                match rng.range(0, 3) {
                    0 => {
                        s.s1 += step;
                        if !allow_unbal { s.p1 += step; }
                    }
                    1 => {
                        s.s2 += step;
                        if !allow_unbal { s.p2 += step; }
                    }
                    2 => {
                        s.p1 += step;
                        if !allow_unbal { s.p2 -= step; }
                    }
                    _ => {
                        s.s1 += step;
                        s.s2 -= step;
                    }
                }
            }
            4 => {
                // Perturb left block if any
                if allow_left && !s.l_blocks.is_empty() {
                    let idx = rng.range(0, s.l_blocks.len() as i32 - 1) as usize;
                    let step = if rng.f01() < 0.5 { 1 } else { -1 };
                    if rng.f01() < 0.5 {
                        s.l_blocks[idx].0 += step;
                    } else {
                        s.l_blocks[idx].1 += step;
                        if !allow_unbal { s.p1 += step; }
                    }
                } else {
                    // Shift two right blocks symmetrically
                    if s.r_blocks.len() >= 2 {
                        let step = if rng.f01() < 0.5 { 1 } else { -1 };
                        let i = rng.range(0, s.r_blocks.len() as i32 - 1) as usize;
                        s.r_blocks[i].0 += step;
                        s.r_blocks[i].1 += step;
                    }
                }
            }
            5 => {
                // Scale-preserving small length transfer between two blocks
                if s.r_blocks.len() >= 2 {
                    let i = rng.range(0, s.r_blocks.len() as i32 - 1) as usize;
                    let j = rng.range(0, s.r_blocks.len() as i32 - 1) as usize;
                    if i != j {
                        s.r_blocks[i].1 += 1;
                        s.r_blocks[j].1 -= 1;
                    }
                }
            }
            6 => {
                // Transfer length between right block and left block!
                if allow_left && !s.l_blocks.is_empty() {
                    let i = rng.range(0, s.r_blocks.len() as i32 - 1) as usize;
                    let j = rng.range(0, s.l_blocks.len() as i32 - 1) as usize;
                    let step = if rng.f01() < 0.5 { 1 } else { -1 };
                    s.r_blocks[i].1 += step;
                    s.l_blocks[j].1 -= step;
                }
            }
            _ => {
                // Global shift s1, s2, p1, p2 (and q_v) by +/-1
                let step = if rng.f01() < 0.5 { 1 } else { -1 };
                if s.s1 + step >= 0 {
                    s.s1 += step; s.s2 += step; s.p1 += step; s.p2 += step;
                    for lb in &mut s.l_blocks { lb.0 += step; }
                }
            }
        }
    }
    // Keep l_blocks q >= p1
    for lb in &mut s.l_blocks {
        if lb.0 < s.p1 { lb.0 = s.p1; }
    }
    s.canonical()
}

fn run_family(name: &str, cands: Vec<GenSpec>, primes: &Arc<Vec<u64>>, top_k: usize) -> Vec<(f64, f64, GenSpec, Saddle)> {
    let t0 = Instant::now();
    // Deduplicate and filter valid saddles
    let mut seen: HashSet<GenSpec> = HashSet::new();
    let mut valid: Vec<(GenSpec, Saddle)> = Vec::new();
    for raw in cands {
        let sp = raw.canonical();
        if !seen.insert(sp.clone()) { continue; }
        if let Some(sr) = eval_saddle(&sp) {
            if sr.m_max <= 48 && sr.g - sr.sigma > 0.45 * (sr.m_max as f64) {
                valid.push((sp, sr));
            }
        }
    }
    let n_valid = valid.len();
    let v_arc = Arc::new(valid);
    let hits_m = Arc::new(Mutex::new(Vec::<(f64, f64, GenSpec, Saddle)>::new()));
    let n_threads = 16;
    let chunk = (n_valid + n_threads - 1) / n_threads;
    let mut handles = Vec::new();
    for tid in 0..n_threads {
        let vr = Arc::clone(&v_arc); let pr = Arc::clone(primes); let hr = Arc::clone(&hits_m);
        handles.push(thread::spawn(move || {
            let lo = tid * chunk; let hi = (lo + chunk).min(vr.len());
            let mut local = Vec::new();
            for i in lo..hi {
                let (ref sp, ref sr) = vr[i];
                if let Some((mu, dt)) = fast_mu(sp, sr, &pr) {
                    if mu < 6.0 { local.push((mu, dt, sp.clone(), sr.clone())); }
                }
            }
            hr.lock().unwrap().extend(local);
        }));
    }
    for h in handles { h.join().unwrap(); }
    let mut hits = Arc::try_unwrap(hits_m).unwrap().into_inner().unwrap();
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("\n[{}] Screened {} valid distinct specs in {:.2}s -> {} with fast_mu < 6.0",
             name, n_valid, t0.elapsed().as_secs_f64(), hits.len());
    // Certify top_k with exact Farey (3 primes/cell)
    let mut certs = Vec::new();
    for (_, _, sp, sr) in hits.iter().take(top_k) {
        if let Some((mu, dt, d1, d2, dis)) = exact_farey_mu(sp, sr, primes, 3) {
            certs.push((mu, dt, d1, d2, dis, sp.clone(), sr.clone()));
        }
    }
    certs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for (rk, (mu, dt, d1, d2, dis, sp, sr)) in certs.iter().take(8).enumerate() {
        let h = sr.sigma + (3 * sr.m_max) as f64 - dt;
        println!("  #{:02} mu={:.6} dis={} def={} | g={:.3} sig={:.3} 3M={} D={:.3}({:.2}+{:.2}) h={:.3} | R={:?} L={:?} P=({},{};{},{})",
                 rk + 1, mu, dis, sp.deg_def(), sr.g, sr.sigma, 3*sr.m_max, dt, d1, d2, h,
                 sp.r_blocks, sp.l_blocks, sp.s1, sp.p1, sp.s2, sp.p2);
    }
    hits.into_iter().take(top_k).collect()
}

fn main() {
    let t_start = Instant::now();
    let primes = Arc::new(sieve_primes(3000));

    // Self-check benchmark: our record F_n*(t) must certify at mu = 5.303682!
    let fstar = GenSpec {
        r_blocks: vec![(0, 13), (2, 15), (4, 17)],
        l_blocks: vec![],
        s1: 2, p1: 23, s2: 3, p2: 21,
    };
    let sr_star = eval_saddle(&fstar).unwrap();
    let (mu_s, dt_s, d1_s, d2_s, dis_s) = exact_farey_mu(&fstar, &sr_star, &primes, 3).unwrap();
    println!("=== SELF-CHECK F_n*(t): mu={:.6} (expected 5.303682), g={:.6}, sig={:.6}, D={:.6} ({:.4}+{:.4}), dis={} ===",
             mu_s, sr_star.g, sr_star.sigma, dt_s, d1_s, d2_s, dis_s);

    // -------------------------------------------------------------------------
    // EXP 6: Coincident zero starts (z1 = z2 = 0 => t0 = 0!) and coincident ends
    // -------------------------------------------------------------------------
    let mut cands_exp6 = Vec::new();
    // 6A: Equal-L coincident starts (0, L), (0, L), (z3, z3+L) for L in 5..=20
    for l in 5..=20 {
        let tot = 3 * l;
        for z3 in 0..=8 {
            for d1 in (tot / 2)..=(tot / 2 + 6).min(tot - 3) {
                let d2 = tot - d1;
                let ds = d1 - d2;
                for s1 in 0..=4 {
                    let s2c = s1 + ds / 2;
                    for s2 in (s2c - 1).max(s1)..=(s2c + 1).min(s1 + ds) {
                        cands_exp6.push(GenSpec {
                            r_blocks: vec![(0, l), (0, l), (z3, z3 + l)],
                            l_blocks: vec![],
                            s1, p1: s1 + d1, s2, p2: s2 + d2,
                        });
                    }
                }
            }
        }
    }
    // 6B: Unequal-L coincident starts (stride 2 on L1, L2, L3)
    for l1 in (6..=16).step_by(2) {
        for l2 in (l1..=16).step_by(2) {
            for l3 in (6..=16).step_by(2) {
                let tot = l1 + l2 + l3;
                for z3 in 0..=5 {
                    for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 3) {
                        let d2 = tot - d1; let ds = d1 - d2;
                        for s1 in 0..=3 {
                            let s2 = s1 + ds / 2;
                            cands_exp6.push(GenSpec {
                                r_blocks: vec![(0, l1), (0, l2), (z3, z3 + l3)],
                                l_blocks: vec![],
                                s1, p1: s1 + d1, s2, p2: s2 + d2,
                            });
                        }
                    }
                }
            }
        }
    }
    // 6C: 4-block with coincident starts (0,0,z3,z4)
    for l in 5..=11 {
        let tot = 4 * l;
        for z3 in 0..=3 {
            for z4 in z3..=(z3 + 3) {
                for d1 in (tot / 2)..=(tot / 2 + 3).min(tot - 3) {
                    let d2 = tot - d1; let ds = d1 - d2;
                    for s1 in 0..=2 {
                        let s2 = s1 + ds / 2;
                        cands_exp6.push(GenSpec {
                            r_blocks: vec![(0, l), (0, l), (z3, z3 + l), (z4, z4 + l)],
                            l_blocks: vec![],
                            s1, p1: s1 + d1, s2, p2: s2 + d2,
                        });
                    }
                }
            }
        }
    }
    run_family("EXP 6: Coincident Zero Starts (t0=0)", cands_exp6, &primes, 15);

    // -------------------------------------------------------------------------
    // EXP 7: Outer pole s1 = 0 (3-block, L in 5..=22)
    // -------------------------------------------------------------------------
    let mut cands_exp7 = Vec::new();
    for l in 5..=22 {
        let tot = 3 * l;
        for d1 in (tot / 2)..=(tot / 2 + 5).min(tot - 3) {
            let d2 = tot - d1; let ds = d1 - d2;
            for s2 in (ds / 2 - 1).max(0)..=((ds + 1) / 2 + 1).min(ds) {
                for z2 in 0..=4 {
                    for z3 in z2..=(z2 + 4) {
                        cands_exp7.push(GenSpec {
                            r_blocks: vec![(0, l), (z2, z2 + l), (z3, z3 + l)],
                            l_blocks: vec![],
                            s1: 0, p1: d1, s2, p2: s2 + d2,
                        });
                    }
                }
            }
        }
    }
    run_family("EXP 7: Outer Pole s1 = 0", cands_exp7, &primes, 15);

    // -------------------------------------------------------------------------
    // EXP 8: Unbalanced / Super-Convergent Degree Excess (Delta_def = 1..8) & 2-Block
    // -------------------------------------------------------------------------
    let mut cands_exp8 = Vec::new();
    for l in 6..=16 {
        for def in 1..=6 {
            let tot = 3 * l + def;
            for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 3) {
                let d2 = tot - d1; let ds = d1 - d2;
                for s1 in 0..=3 {
                    let s2 = s1 + ds / 2;
                    for z2 in 0..=3 {
                        for z3 in z2..=(z2 + 3) {
                            cands_exp8.push(GenSpec {
                                r_blocks: vec![(0, l), (z2, z2 + l), (z3, z3 + l)],
                                l_blocks: vec![],
                                s1, p1: s1 + d1, s2, p2: s2 + d2,
                            });
                        }
                    }
                }
            }
        }
    }
    // Also 2-block numerator (m_R = 2, balanced & unbalanced!)
    for l1 in 5..=16 {
        for l2 in l1..=16 {
            for def in 0..=8 {
                let tot = l1 + l2 + def;
                for d1 in (tot / 2)..=(tot / 2 + 4).min(tot - 2) {
                    let d2 = tot - d1; let ds = d1 - d2;
                    for s1 in 0..=3 {
                        let s2 = s1 + ds / 2;
                        for z2 in 0..=3 {
                            cands_exp8.push(GenSpec {
                                r_blocks: vec![(0, l1), (z2, z2 + l2)],
                                l_blocks: vec![],
                                s1, p1: s1 + d1, s2, p2: s2 + d2,
                            });
                        }
                    }
                }
            }
        }
    }
    run_family("EXP 8: Unbalanced Super-Convergent (Delta_def > 0) & 2-Block", cands_exp8, &primes, 15);

    // -------------------------------------------------------------------------
    // EXP 21: TWO-SIDED ZERO BLOCKS (Right zeros + LEFT zeros at t < -p1*n!)
    // -------------------------------------------------------------------------
    let mut cands_exp21 = Vec::new();
    // 21A: m_R = 2, m_L = 1 (still 3 total blocks -> 5-binomial 5F4!)
    for lr in 5..=15 {
        for kl in 3..=15 {
            for def in [0, 1] {
                let tot = 2 * lr + kl + def;
                if tot > 42 { continue; }
                for d1 in (tot / 2)..=(tot / 2 + 3).min(tot - 3) {
                    let d2 = tot - d1; let ds = d1 - d2;
                    for s1 in 0..=3 {
                        let s2 = s1 + ds / 2;
                        let p1 = s1 + d1; let p2 = s2 + d2;
                        for z2 in 0..=3 {
                            for dq in 0..=3 {
                                cands_exp21.push(GenSpec {
                                    r_blocks: vec![(0, lr), (z2, z2 + lr)],
                                    l_blocks: vec![(p1 + dq, kl)],
                                    s1, p1, s2, p2,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    // 21B: m_R = 2, m_L = 2 (symmetric 4-block two-sided!)
    for lr in 5..=11 {
        for kl in 4..=11 {
            let tot = 2 * lr + 2 * kl;
            if tot > 42 { continue; }
            for d1 in (tot / 2)..=(tot / 2 + 3).min(tot - 3) {
                let d2 = tot - d1; let ds = d1 - d2;
                for s1 in 0..=2 {
                    let s2 = s1 + ds / 2;
                    let p1 = s1 + d1; let p2 = s2 + d2;
                    for z2 in 0..=2 {
                        for dq1 in 0..=2 {
                            for dq2 in dq1..=(dq1 + 2) {
                                cands_exp21.push(GenSpec {
                                    r_blocks: vec![(0, lr), (z2, z2 + lr)],
                                    l_blocks: vec![(p1 + dq1, kl), (p1 + dq2, kl)],
                                    s1, p1, s2, p2,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    // 21C: m_R = 3, m_L = 1 (3-block right + 1 left stabilizer!)
    for lr in 6..=12 {
        for kl in 2..=8 {
            let tot = 3 * lr + kl;
            if tot > 42 { continue; }
            for d1 in (tot / 2)..=(tot / 2 + 3).min(tot - 3) {
                let d2 = tot - d1; let ds = d1 - d2;
                for s1 in 1..=2 {
                    let s2 = s1 + ds / 2;
                    let p1 = s1 + d1; let p2 = s2 + d2;
                    for (z2, z3) in [(1, 2), (2, 4), (0, 2)] {
                        for dq in 0..=2 {
                            cands_exp21.push(GenSpec {
                                r_blocks: vec![(0, lr), (z2, z2 + lr), (z3, z3 + lr)],
                                l_blocks: vec![(p1 + dq, kl)],
                                s1, p1, s2, p2,
                            });
                        }
                    }
                }
            }
        }
    }
    let best_21 = run_family("EXP 21: TWO-SIDED ZERO BLOCKS (Right + Left Zeros!)", cands_exp21, &primes, 20);

    // -------------------------------------------------------------------------
    // EXP 9 & EXP 10: 16-Core Multi-Start Evolutionary / Simulated Annealing
    // Optimizer across m_R in {3, 4, 5, 6}, m_L in {0, 1, 2}, high-L & unequal-L!
    // -------------------------------------------------------------------------
    println!("\n=== EXP 9 & 10: 16-CORE MULTI-START EVOLUTIONARY OPTIMIZER (m=3,4,5,6; 1-sided & 2-sided; balanced & unbalanced) ===");
    let t_evo = Instant::now();
    // Build diverse seed pool across all topologies!
    let mut seeds: Vec<(String, GenSpec, bool, bool)> = vec![
        ("3blk-Fstar".into(), fstar.clone(), false, false),
        ("3blk-Fstar-unbal".into(), fstar.clone(), true, false),
        ("3blk-2xFstar".into(), GenSpec { r_blocks: vec![(0,26),(4,30),(8,34)], l_blocks: vec![], s1: 4, p1: 46, s2: 6, p2: 42 }, false, false),
        ("3blk-LK".into(), GenSpec { r_blocks: vec![(0,10),(1,11),(2,12)], l_blocks: vec![], s1: 1, p1: 17, s2: 2, p2: 16 }, false, false),
        ("3blk-highL".into(), GenSpec { r_blocks: vec![(0,18),(3,21),(6,24)], l_blocks: vec![], s1: 3, p1: 32, s2: 5, p2: 30 }, false, false),
        ("3blk-wide".into(), GenSpec { r_blocks: vec![(0,15),(2,19),(5,22)], l_blocks: vec![], s1: 2, p1: 28, s2: 4, p2: 26 }, true, false),
        ("4blk-eq".into(), GenSpec { r_blocks: vec![(0,10),(1,11),(2,12),(3,13)], l_blocks: vec![], s1: 2, p1: 24, s2: 4, p2: 22 }, false, false),
        ("4blk-unbal".into(), GenSpec { r_blocks: vec![(0,9),(2,11),(4,13),(6,15)], l_blocks: vec![], s1: 2, p1: 22, s2: 3, p2: 21 }, true, false),
        ("4blk-coinc".into(), GenSpec { r_blocks: vec![(0,11),(0,13),(2,15),(4,17)], l_blocks: vec![], s1: 1, p1: 27, s2: 3, p2: 26 }, false, false),
        ("5blk-eq".into(), GenSpec { r_blocks: vec![(0,8),(1,9),(2,10),(3,11),(4,12)], l_blocks: vec![], s1: 2, p1: 23, s2: 3, p2: 22 }, false, false),
        ("5blk-unbal".into(), GenSpec { r_blocks: vec![(0,7),(0,8),(2,9),(3,10),(5,12)], l_blocks: vec![], s1: 1, p1: 21, s2: 2, p2: 20 }, true, false),
        ("6blk-eq".into(), GenSpec { r_blocks: vec![(0,6),(1,7),(2,8),(3,9),(4,10),(5,11)], l_blocks: vec![], s1: 1, p1: 20, s2: 2, p2: 19 }, false, false),
        ("2side-2R1L".into(), GenSpec { r_blocks: vec![(0,12),(2,14)], l_blocks: vec![(21,10)], s1: 2, p1: 20, s2: 3, p2: 19 }, false, true),
        ("2side-2R2L".into(), GenSpec { r_blocks: vec![(0,10),(2,12)], l_blocks: vec![(22,10),(24,10)], s1: 2, p1: 22, s2: 3, p2: 21 }, false, true),
        ("2side-3R1L".into(), GenSpec { r_blocks: vec![(0,11),(2,13),(4,15)], l_blocks: vec![(23,6)], s1: 2, p1: 22, s2: 3, p2: 21 }, false, true),
        ("2side-3R2L".into(), GenSpec { r_blocks: vec![(0,9),(2,11),(4,13)], l_blocks: vec![(22,7),(24,7)], s1: 2, p1: 22, s2: 3, p2: 21 }, true, true),
    ];
    // Replace 2-sided seeds with best_21 if any beat them!
    if let Some((_, _, sp0, _)) = best_21.first() {
        seeds[12].1 = sp0.clone();
    }

    let evo_res = Arc::new(Mutex::new(Vec::<(f64, String, GenSpec, Saddle)>::new()));
    let mut handles = Vec::new();
    for (tid, (label, seed, allow_unbal, allow_left)) in seeds.into_iter().enumerate() {
        let pr = Arc::clone(&primes);
        let er = Arc::clone(&evo_res);
        handles.push(thread::spawn(move || {
            let mut rng = Rng(0x9E3779B97F4A7C15u64 ^ ((tid as u64 + 1) * 0xBF58476D1CE4E5B9));
            let mut cur = seed.canonical();
            let mut cur_sr = eval_saddle(&cur).unwrap_or_else(|| eval_saddle(&GenSpec {
                r_blocks: vec![(0,10),(1,11),(2,12)], l_blocks: vec![], s1: 1, p1: 17, s2: 2, p2: 16
            }).unwrap());
            let mut cur_mu = fast_mu(&cur, &cur_sr, &pr).map(|x| x.0).unwrap_or(10.0);
            let mut best = (cur_mu, cur.clone(), cur_sr.clone());
            let mut top_local: Vec<(f64, GenSpec, Saddle)> = vec![best.clone()];
            let mut seen_loc: HashSet<GenSpec> = HashSet::new();
            seen_loc.insert(cur.clone());

            let steps = 1600;
            for st in 0..steps {
                let temp = 0.25 * (0.015_f64).powf(st as f64 / steps as f64);
                let cand = mutate_spec(&cur, &mut rng, allow_unbal, allow_left);
                if let Some(sr) = eval_saddle(&cand) {
                    if sr.m_max > 56 || max_arg_1(&cand, &sr) > 110 { continue; }
                    if let Some((mu, _)) = fast_mu(&cand, &sr, &pr) {
                        let accept = mu < cur_mu || rng.f01() < ((cur_mu - mu) / temp).exp();
                        if accept {
                            cur = cand.clone(); cur_sr = sr.clone(); cur_mu = mu;
                        }
                        if mu < 5.65 && seen_loc.insert(cand.clone()) {
                            top_local.push((mu, cand.clone(), sr.clone()));
                        }
                        if mu < best.0 {
                            best = (mu, cand, sr);
                        }
                    }
                }
            }
            top_local.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let mut lock = er.lock().unwrap();
            for (mu, sp, sr) in top_local.into_iter().take(6) {
                lock.push((mu, label.clone(), sp, sr));
            }
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut evo_hits = Arc::try_unwrap(evo_res).unwrap().into_inner().unwrap();
    evo_hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut seen_evo = HashSet::new();
    let mut unique_evo = Vec::new();
    for (mu, lbl, sp, sr) in evo_hits {
        if seen_evo.insert(sp.clone()) { unique_evo.push((mu, lbl, sp, sr)); }
    }
    println!("Evolutionary search finished in {:.2}s. Certifying top 25 unique winners with EXACT Farey (3 primes/cell)...",
             t_evo.elapsed().as_secs_f64());
    let mut final_certs = Vec::new();
    for (_, lbl, sp, sr) in unique_evo.into_iter().take(25) {
        if let Some((mu, dt, d1, d2, dis)) = exact_farey_mu(&sp, &sr, &primes, 3) {
            final_certs.push((mu, dt, d1, d2, dis, lbl, sp, sr));
        }
    }
    final_certs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("\n=== TOP 20 EVOLUTIONARY WINNERS (EXACT FAREY CERTIFIED) ===");
    for (rk, (mu, dt, _d1, _d2, dis, lbl, sp, sr)) in final_certs.iter().take(20).enumerate() {
        let h = sr.sigma + (3 * sr.m_max) as f64 - dt;
        println!("  #{:02} [{:16}] mu={:.6} dis={} def={} | g={:.3} sig={:.3} 3M={} D={:.3} h={:.3} | R={:?} L={:?} P=({},{};{},{})",
                 rk + 1, lbl, mu, dis, sp.deg_def(), sr.g, sr.sigma, 3*sr.m_max, dt, h,
                 sp.r_blocks, sp.l_blocks, sp.s1, sp.p1, sp.s2, sp.p2);
    }
    println!("\nTotal Phase 2 elapsed: {:.2}s.", t_start.elapsed().as_secs_f64());
}
