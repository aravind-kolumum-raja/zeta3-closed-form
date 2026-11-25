/-!
# Auxiliary Lean 4 Algebraic and Arithmetic Checks for F_n*(t)

Zero-dependency Lean 4 (v4.36.0) companion checks for:

  "Improved Irrationality Measure mu(zeta(3)) < 5.3037 and Single-Presentation Linear
   Independence of 1, zeta(3), zeta(4) via Elementary Congruences"

Scope: this file machine-checks auxiliary algebraic identities (the saddle cubic,
1-periodicity of Delta, the Wilson-Lucas degree identity, smooth Z/p^2Z crossing
identities, sum(beta) = -1, quotient degree balance of G_tilde_q), the strict
concavity polynomial at c = 43/4, the 38 first-order and 5 second-order sub-shell
midpoint defects, the exact rational u < 1 sums, and the final cross-multiplied
integer inequalities. The analytic, asymptotic, congruence, and Prime Number Theorem
proofs themselves are carried out in the paper and tested by standalone Python scripts.

Zero `sorry`, zero `native_decide`, zero extra axioms (pure Lean 4 kernel).
-/

namespace Zeta3.ElementaryRecordProof

/-! ## 1. Exact Primitive Saddle Cubic of F_n*(t) -/

theorem saddle_cubic_identity (x : Int) :
    x * (x - 2) * (x - 4) * (x + 2) * (x + 3)
      - (x - 13) * (x - 15) * (x - 17) * (x + 23) * (x + 21)
    = 15 * (54 * x ^ 3 - 298 * x ^ 2 - 11879 * x + 106743) := by
  grind

def saddleCubic (x : Int) : Int :=
  54 * x ^ 3 - 298 * x ^ 2 - 11879 * x + 106743

theorem saddle_pole_root_bracket :
    And (saddleCubic (-16) < 0) (saddleCubic (-15) > 0) := by
  decide

/-! ## 2. Exact Wilson-Lucas Degree Identity & Periodicity of Delta(u) -/

def ratFloor (a : Int) (b : Nat) : Int :=
  a / (b : Int)

def wilsonLucasDelta (a : Int) (b : Nat) : Int :=
  2 - 3 * ratFloor (13 * a) b + ratFloor (20 * a) b + ratFloor (19 * a) b

theorem wilsonLucasDelta_periodic (a : Int) (b : Nat) (hb : b > 0) :
    wilsonLucasDelta (a + b) b = wilsonLucasDelta a b := by
  unfold wilsonLucasDelta ratFloor
  have hbp : Ne (b : Int) 0 := by omega
  have h13 : 13 * (a + (b : Int)) = 13 * a + (b : Int) * 13 := by omega
  have h20 : 20 * (a + (b : Int)) = 20 * a + (b : Int) * 20 := by omega
  have h19 : 19 * (a + (b : Int)) = 19 * a + (b : Int) * 19 := by omega
  rw [h13, h20, h19,
      Int.add_mul_ediv_left _ _ hbp,
      Int.add_mul_ediv_left _ _ hbp,
      Int.add_mul_ediv_left _ _ hbp]
  omega

theorem wilson_lucas_degree_eq_p_minus_2
    (p n q13 r13 q20 r20 q19 r19 c1 c2 c3 c4 c5 : Int)
    (h13 : 13 * n = q13 * p + r13)
    (h20 : 20 * n = q20 * p + r20)
    (h19 : 19 * n = q19 * p + r19)
    (h_defect : (2 - 3 * q13 + q20 + q19) - (c1 + c2 + c3 + c4 + c5) = 1) :
    (r13 - c1 * p) + (r13 - c2 * p) + (r13 - c3 * p)
      + ((p - 1 - r20) - c4 * p) + ((p - 1 - r19) - c5 * p) = p - 2 := by
  grind

/-! ## 3. Exact Midpoint Defect Check for All 38 Uniform-Carry Sub-Shells -/

structure SubShell where
  an : Nat
  ad : Nat
  bn : Nat
  bd : Nat
  vmin : Int

def subShells38 : List SubShell :=
  [ SubShell.mk 1 15 2 27 3,   SubShell.mk 1 13 2 25 0,   SubShell.mk 1 11 2 21 0,
    SubShell.mk 2 17 3 25 2,   SubShell.mk 2 13 3 19 0,   SubShell.mk 3 19 4 25 1,
    SubShell.mk 3 16 4 21 1,   SubShell.mk 3 13 4 17 0,   SubShell.mk 4 17 5 21 0,
    SubShell.mk 5 21 6 25 0,   SubShell.mk 3 11 7 25 2,   SubShell.mk 5 16 6 19 0,
    SubShell.mk 6 19 8 25 1,   SubShell.mk 5 13 9 23 0,   SubShell.mk 6 13 7 15 0,
    SubShell.mk 7 15 8 17 0,   SubShell.mk 8 17 9 19 0,   SubShell.mk 9 19 10 21 1,
    SubShell.mk 10 21 11 23 1, SubShell.mk 11 23 12 25 1, SubShell.mk 7 13 11 20 0,
    SubShell.mk 11 20 14 25 1, SubShell.mk 14 25 13 23 1, SubShell.mk 8 13 13 21 0,
    SubShell.mk 11 17 13 20 1, SubShell.mk 13 20 15 23 2, SubShell.mk 9 13 16 23 0,
    SubShell.mk 11 15 14 19 1, SubShell.mk 14 19 17 23 2, SubShell.mk 10 13 18 23 0,
    SubShell.mk 14 17 19 23 2, SubShell.mk 11 13 17 20 0, SubShell.mk 17 20 23 27 1,
    SubShell.mk 14 15 16 17 0, SubShell.mk 16 17 18 19 0, SubShell.mk 18 19 19 20 1,
    SubShell.mk 19 20 20 21 2, SubShell.mk 20 21 22 23 2 ]

def checkSubShell (s : SubShell) : Bool :=
  let ordered := s.an * s.bd < s.bn * s.ad
  let midN : Int := (s.an * s.bd + s.bn * s.ad : Nat)
  let midD : Nat := 2 * s.ad * s.bd
  ordered && (wilsonLucasDelta midN midD - s.vmin == 1)

theorem all_38_subshells_defect_eq_one :
    And (subShells38.length = 38) (subShells38.all checkSubShell = true) := by
  decide

/-! ## 4. Algebraic Strict Concavity Polynomial of psi''(y) at c = 43/4 -/

def polyMul (p q : List Int) : List Int :=
  let n := p.length + q.length - 1
  List.ofFn (fun (idx : Fin n) =>
    let k := idx.val
    (List.range p.length).foldl (fun acc i =>
      if And (i <= k) (k - i < q.length) then
        acc + (p.getD i 0) * (q.getD (k - i) 0)
      else acc) 0)

def contourTermsAt43_4 : List (Prod Int Int) :=
  [ (-1, 43), (1, -9), (-1, 35), (1, -17), (-1, 27), (1, -25),
    (1, 135), (-1, 51), (1, 127), (-1, 55) ]

def contourConcavityCoeffs : List Int :=
  let n := contourTermsAt43_4.length
  let init := List.replicate n (0 : Int)
  (List.range n).foldl (fun acc i =>
    let (eps, ai) := contourTermsAt43_4.getD i (0, 0)
    let prod := (List.range n).foldl (fun pr j =>
      if j == i then pr
      else
        let (_, aj) := contourTermsAt43_4.getD j (0, 0)
        polyMul pr [aj * aj, 1]) [1]
    List.zipWith (fun x y => x + eps * ai * y) acc prod) init

theorem contour_strictly_concave_at_43_4 :
    And (contourConcavityCoeffs.getD 9 1 = 0)
        ((contourConcavityCoeffs.take 9).all (fun c => c < 0) = true) := by
  decide

/-! ## 5. Second-Order Z/p^2Z Crossing Identities, sum(beta) = -1, Quotient Balance -/

def binom : Nat -> Nat -> Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => binom n k + binom n (k + 1)

theorem binom_eq_zero_of_lt : forall {n k : Nat}, n < k -> binom n k = 0
  | 0, _ + 1, _ => rfl
  | n + 1, k + 1, h => by
    simp only [binom]
    rw [binom_eq_zero_of_lt (by omega : n < k), binom_eq_zero_of_lt (by omega : n < k + 1)]
  | _, 0, h => by omega

theorem num_crossing_identity : forall (q a : Nat),
    (q + 1 - a) * binom (q + 1) a = (q + 1) * binom q a
  | 0, 0 => rfl
  | 0, 1 => rfl
  | 0, _ + 2 => rfl
  | q + 1, 0 => by simp [binom]
  | q + 1, a' + 1 => by
    by_cases hle : a' <= q
    case pos =>
      have h1 := num_crossing_identity q a'
      have h2 := num_crossing_identity q (a' + 1)
      have hp1 : binom (q + 2) (a' + 1) = binom (q + 1) a' + binom (q + 1) (a' + 1) := rfl
      have hp2 : binom (q + 1) (a' + 1) = binom q a' + binom q (a' + 1) := rfl
      have hsub1 : q + 1 + 1 - (a' + 1) = q + 1 - a' := by omega
      have hsub2 : q + 1 - a' = (q + 1 - (a' + 1)) + 1 := by omega
      rw [hp1, hsub1, Nat.mul_add, h1]
      have hterm : (q + 1 - a') * binom (q + 1) (a' + 1)
          = (q + 1) * binom q (a' + 1) + binom (q + 1) (a' + 1) := by
        rw [hsub2, Nat.add_mul, Nat.one_mul, h2]
      rw [hterm, hp2]
      grind
    case neg =>
      by_cases heq : a' = q + 1
      case pos =>
        subst heq
        have hzero : binom (q + 1) (q + 2) = 0 := binom_eq_zero_of_lt (by omega)
        simp [hzero]
      case neg =>
        have hgt1 : q + 2 < a' + 1 := by omega
        have hgt2 : q + 1 < a' + 1 := by omega
        rw [binom_eq_zero_of_lt hgt1, binom_eq_zero_of_lt hgt2]
        simp

theorem den_crossing_identity : forall (e k : Nat),
    (e - k) * binom e k = (k + 1) * binom e (k + 1)
  | 0, _ => by simp [binom]
  | e + 1, 0 => by
    have h1 := den_crossing_identity e 0
    simp only [binom, Nat.sub_zero, Nat.mul_one] at h1
    simp only [binom, Nat.sub_zero, Nat.mul_one]
    omega
  | e + 1, k' + 1 => by
    by_cases hle : k' + 1 <= e
    case pos =>
      have h1 := den_crossing_identity e k'
      have h2 := den_crossing_identity e (k' + 1)
      have hp1 : binom (e + 1) (k' + 1) = binom e k' + binom e (k' + 1) := rfl
      have hp2 : binom (e + 1) (k' + 2) = binom e (k' + 1) + binom e (k' + 2) := rfl
      have hsub1 : e + 1 - (k' + 1) = e - k' := by omega
      have hsub2 : e - k' = (e - (k' + 1)) + 1 := by omega
      rw [hp1, hp2, hsub1, Nat.mul_add, h1]
      have hterm : (e - k') * binom e (k' + 1)
          = (k' + 2) * binom e (k' + 2) + binom e (k' + 1) := by
        rw [hsub2, Nat.add_mul, Nat.one_mul, h2]
      rw [hterm]
      grind
    case neg =>
      by_cases heq : k' = e
      case pos =>
        subst heq
        have hzero : binom (k' + 1) (k' + 2) = 0 := binom_eq_zero_of_lt (by omega)
        simp [hzero]
      case neg =>
        have hgt1 : e + 1 < k' + 1 := by omega
        have hgt2 : e + 1 < k' + 2 := by omega
        rw [binom_eq_zero_of_lt hgt1, binom_eq_zero_of_lt hgt2]
        simp

theorem second_order_sum_beta_eq_neg_one
    (q13 q20 q19 c4 c5 : Int)
    (h_defect : (2 - 3 * q13 + q20 + q19) - (c4 + c5) = 1) :
    3 * q13 - (q20 + 1 - c4) - (q19 + 1 - c5) = -1 := by
  grind

theorem second_order_quotient_G_tilde_balanced (q : Nat) :
    And ((13 * q + 12) + (13 * q + 12) + (13 * q + 13) = 39 * q + 37)
    (And ((20 * q + 19) + (19 * q + 18) = 39 * q + 37)
    (And (4 * q + 4 >= 1) (2 * q + 2 >= 1))) := by
  omega

def subWindows2ndOrder5 : List SubShell :=
  [ SubShell.mk 7 15 9 19 0,   SubShell.mk 9 19 11 23 1,
    SubShell.mk 16 17 18 19 0, SubShell.mk 18 19 19 20 1, SubShell.mk 19 20 20 21 2 ]

theorem all_5_second_order_subwindows_defect_eq_one :
    And (subWindows2ndOrder5.length = 5)
        (subWindows2ndOrder5.all checkSubShell = true) := by
  decide

/-! ## 6. Exact Rational u < 1 Sums and Final Cross-Multiplied Integer Inequalities -/

theorem exact_u1_carry_joint_minus_z3 :
    (1193959693 : Int) - 1019365273 = 6 * 29099070 := by
  decide

def L23 : Nat := 5354228880

def subShell38PosSum : Nat :=
  subShells38.foldl (fun acc s => acc + (L23 / s.an) * s.ad) 0

def subShell38NegSum : Nat :=
  subShells38.foldl (fun acc s => acc + (L23 / s.bn) * s.bd) 0

theorem subshells38_exact_u1_sum :
    And (subShells38.all (fun s => L23 % s.an == 0 && L23 % s.bn == 0) = true)
        (subShell38PosSum = subShell38NegSum + 19552742437) := by
  decide

theorem second_order_exact_u1_sum :
    And ((15 * 880 + 17 * 385) = (23 * 560 + 21 * 308) + 397)
        (19552742437 + 397 * (5354228880 / 6160) = 2 * 9948906029) := by
  decide

theorem thm1_mu_zeta3_lt_53037 :
    150301021 * 10000 < 53037 * (150301021 - 121962024) := by
  decide

theorem thm1_mu_zeta3_lt_53193 :
    150301021 * 10000 < 53193 * (150301021 - 122045057) := by
  decide

theorem thm2_mu_joint_lt_112_3 :
    And (148962024 < 150301021)
        (150301021 * 10 < 1123 * (150301021 - 148962024)) := by
  decide

theorem thm2_mu_joint_lt_119_7 :
    And (149045057 < 150301021)
        (150301021 * 10 < 1197 * (150301021 - 149045057)) := by
  decide

end Zeta3.ElementaryRecordProof
