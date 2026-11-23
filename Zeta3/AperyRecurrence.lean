import Zeta3.Basic

/-!
# Apéry's Recurrence & The Casoratian Conservation Theorem

We prove completely formally (with zero `sorry` and no external axioms):
1. `aperySeq_satisfies_recurrence`: `aperySeq u0 u1` satisfies Apéry's recurrence.
2. `casoratian_step`: For ANY two solutions `u, v` of Apéry's recurrence,
   $$(m+2)^3 (u_{m+1} v_{m+2} - u_{m+2} v_{m+1}) = (m+1)^3 (u_m v_{m+1} - u_{m+1} v_m).$$
3. `casoratian_const`: By induction, for all `n : Nat`,
   $$(n+1)^3 (u_n v_{n+1} - u_{n+1} v_n) = u_0 v_1 - u_1 v_0.$$
4. `apery_casoratian`: Specializing to `aperyA` and `aperyB`,
   $$(n+1)^3 (a_n b_{n+1} - a_{n+1} b_n) = 6 \qquad (\forall n \in \mathbb{N}).$$
-/

namespace Zeta3

/-- Helper: `(k + 1 : Rat) ≠ 0`. -/
theorem nat_succ_rat_ne_zero (k : Nat) : (((k + 1 : Nat) : Rat)) ≠ 0 := by
  grind

/-- Helper: `((k + 1 : Nat) : Rat) ^ 3 ≠ 0`. -/
theorem cube_succ_rat_ne_zero (k : Nat) : (((k + 1 : Nat) : Rat) ^ 3) ≠ 0 := by
  have h := nat_succ_rat_ne_zero k
  intro hcube
  have hmul :
      (((k + 1 : Nat) : Rat) * ((k + 1 : Nat) : Rat) * ((k + 1 : Nat) : Rat)) = 0 := by
    grind
  rcases Rat.mul_eq_zero.mp hmul with h12 | h3
  · rcases Rat.mul_eq_zero.mp h12 with h1 | h2
    · exact h h1
    · exact h h2
  · exact h h3

/-- Algebraic identity for eliminating the middle term `P * u_{m+1} * v_{m+1}` in a 3-term recurrence. -/
theorem casoratian_alg_identity (c2 c0 P u0 u1 u2 v0 v1 v2 : Rat)
    (hu : c2 * u2 - P * u1 + c0 * u0 = 0)
    (hv : c2 * v2 - P * v1 + c0 * v0 = 0) :
    c2 * (u1 * v2 - u2 * v1) = c0 * (u0 * v1 - u1 * v0) := by
  grind

/--
**Casoratian Step Lemma:**
For any two solutions `u, v` of Apéry's recurrence, the weighted Casoratian
`(m + 1)^3 * (u m * v (m + 1) - u (m + 1) * v m)` is invariant under `m ↦ m + 1`.
-/
theorem casoratian_step (u v : Nat → Rat)
    (hu : SatisfiesAperyRecurrence u) (hv : SatisfiesAperyRecurrence v) (m : Nat) :
    (((m + 2 : Nat) : Rat) ^ 3) * (u (m + 1) * v (m + 2) - u (m + 2) * v (m + 1)) =
    (((m + 1 : Nat) : Rat) ^ 3) * (u m * v (m + 1) - u (m + 1) * v m) :=
  casoratian_alg_identity
    (((m + 2 : Nat) : Rat) ^ 3)
    (((m + 1 : Nat) : Rat) ^ 3)
    (aperyP (m + 1) : Rat)
    (u m) (u (m + 1)) (u (m + 2))
    (v m) (v (m + 1)) (v (m + 2))
    (hu m) (hv m)

/--
**Casoratian Conservation Theorem:**
For any two solutions `u, v` of Apéry's recurrence and every `n : Nat`,
$$(n + 1)^3 (u_n v_{n+1} - u_{n+1} v_n) = u_0 v_1 - u_1 v_0.$$
-/
theorem casoratian_const (u v : Nat → Rat)
    (hu : SatisfiesAperyRecurrence u) (hv : SatisfiesAperyRecurrence v) :
    ∀ n : Nat,
      (((n + 1 : Nat) : Rat) ^ 3) * (u n * v (n + 1) - u (n + 1) * v n) =
      u 0 * v 1 - u 1 * v 0 := by
  intro n
  induction n with
  | zero =>
      grind
  | succ m ih =>
      rw [casoratian_step u v hu hv m, ih]

/-- `aperySeq u0 u1` satisfies Apéry's recurrence. -/
theorem aperySeq_satisfies_recurrence (u0 u1 : Rat) :
    SatisfiesAperyRecurrence (aperySeq u0 u1) := by
  intro m
  simp only [aperySeq]
  have hc2 : (((m + 2 : Nat) : Rat) ^ 3) ≠ 0 := cube_succ_rat_ne_zero (m + 1)
  have hcancel :
      ∀ (a b : Rat), b ≠ 0 → b * (a / b) = a := by
    intro a b hb
    rw [Rat.div_def, ← Rat.mul_assoc, Rat.mul_comm b a, Rat.mul_assoc,
        Rat.mul_inv_cancel b hb, Rat.mul_one]
  rw [hcancel _ _ hc2]
  grind

/--
**Apéry's Casoratian Identity (1978):**
For all `n : Nat`,
$$(n + 1)^3 (a_n b_{n+1} - a_{n+1} b_n) = 6.$$
-/
theorem apery_casoratian (n : Nat) :
    (((n + 1 : Nat) : Rat) ^ 3) * (aperyA n * aperyB (n + 1) - aperyA (n + 1) * aperyB n) = 6 := by
  have h := casoratian_const aperyA aperyB
    (aperySeq_satisfies_recurrence 1 5)
    (aperySeq_satisfies_recurrence 0 6) n
  have h01 : aperyA 0 * aperyB 1 - aperyA 1 * aperyB 0 = 6 := by
    simp only [aperyA, aperyB, aperySeq]
    decide +kernel
  rw [h, h01]

/--
**Apéry Ratio Step Identity:**
Whenever `aperyA n ≠ 0` and `aperyA (n + 1) ≠ 0`,
$$\frac{b_{n+1}}{a_{n+1}} - \frac{b_n}{a_n} = \frac{6}{(n+1)^3 a_n a_{n+1}}.$$
-/
theorem apery_ratio_step (n : Nat) (hn : aperyA n ≠ 0) (hn1 : aperyA (n + 1) ≠ 0) :
    aperyB (n + 1) / aperyA (n + 1) - aperyB n / aperyA n =
    6 / ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) := by
  have hcas := apery_casoratian n
  have hcube : (((n + 1 : Nat) : Rat) ^ 3) ≠ 0 := cube_succ_rat_ne_zero n
  have hprod : (((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1) ≠ 0 := by
    intro h0
    rcases Rat.mul_eq_zero.mp h0 with hleft | hright
    · rcases Rat.mul_eq_zero.mp hleft with hc | ha
      · exact hcube hc
      · exact hn ha
    · exact hn1 hright
  -- Clear denominators using field arithmetic on Rat
  have hmul :
      (aperyB (n + 1) / aperyA (n + 1) - aperyB n / aperyA n) *
        ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) = 6 := by
    have hdiv_cancel : ∀ (x y : Rat), y ≠ 0 → (x / y) * y = x := by
      intro x y hy
      rw [Rat.div_def, Rat.mul_assoc, Rat.inv_mul_cancel y hy, Rat.mul_one]
    have h1 : (aperyB (n + 1) / aperyA (n + 1)) *
        ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) =
        (((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyB (n + 1) := by
      have hrw : (aperyB (n + 1) / aperyA (n + 1)) *
          ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) =
          (((n + 1 : Nat) : Rat) ^ 3) * aperyA n * ((aperyB (n + 1) / aperyA (n + 1)) * aperyA (n + 1)) := by
        grind
      rw [hrw, hdiv_cancel _ _ hn1]
    have h2 : (aperyB n / aperyA n) *
        ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) =
        (((n + 1 : Nat) : Rat) ^ 3) * aperyA (n + 1) * aperyB n := by
      have hrw : (aperyB n / aperyA n) *
          ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) =
          (((n + 1 : Nat) : Rat) ^ 3) * aperyA (n + 1) * ((aperyB n / aperyA n) * aperyA n) := by
        grind
      rw [hrw, hdiv_cancel _ _ hn]
    have hsub :
        (aperyB (n + 1) / aperyA (n + 1) - aperyB n / aperyA n) *
          ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) =
        (aperyB (n + 1) / aperyA (n + 1)) *
          ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) -
        (aperyB n / aperyA n) *
          ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1)) := by
      grind
    rw [hsub, h1, h2]
    grind
  have hdiv : ∀ (x y z : Rat), y ≠ 0 → x * y = z → x = z / y := by
    intro x y z hy hxy
    rw [← hxy, Rat.div_def, Rat.mul_assoc, Rat.mul_inv_cancel y hy, Rat.mul_one]
  exact hdiv _ _ _ hprod hmul

/-- Helper: finite telescoping sum `∑_{k=0}^{N-1} 6 / ((k+1)^3 * a_k * a_{k+1})`. -/
def aperyTelescopingSum : Nat → Rat
  | 0 => 0
  | n + 1 =>
      aperyTelescopingSum n +
        6 / ((((n + 1 : Nat) : Rat) ^ 3) * aperyA n * aperyA (n + 1))

/--
**Apéry's Telescoping Series Theorem:**
Assuming `aperyA k ≠ 0` for all `k ≤ N`, the `N`-th Apéry rational approximant equals
$$\frac{b_N}{a_N} = \sum_{k=0}^{N-1} \frac{6}{(k+1)^3 a_k a_{k+1}}.$$
(Taking `N → ∞` yields Apéry's celebrated series $\zeta(3) = \sum_{k=0}^\infty \frac{6}{(k+1)^3 a_k a_{k+1}}$.)
-/
theorem apery_telescoping_identity :
    ∀ N : Nat, (∀ k ≤ N, aperyA k ≠ 0) →
      aperyB N / aperyA N = aperyTelescopingSum N := by
  intro N
  induction N with
  | zero =>
      intro _
      simp only [aperyTelescopingSum, aperyA, aperyB, aperySeq]
      decide +kernel
  | succ m ih =>
      intro hne
      have ihm := ih (fun k hk => hne k (Nat.le_succ_of_le hk))
      have hstep := apery_ratio_step m (hne m (Nat.le_succ m)) (hne (m + 1) (Nat.le_refl _))
      simp only [aperyTelescopingSum, ← ihm, ← hstep]
      grind

end Zeta3
