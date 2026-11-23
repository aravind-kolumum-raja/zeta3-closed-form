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

end Zeta3
