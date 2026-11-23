/-!
# Euler's Parity Obstruction for Odd Zeta Values

In contour integration / Mittag-Leffler expansions of $\pi \cot(\pi z) / z^m$,
residues at the symmetric integer poles $+n$ and $-n$ pair up as
$\frac{1}{n^m} + \frac{1}{(-n)^m}$.

This module formally proves (in pure Lean 4, zero `sorry`, standard axioms only):
* `neg_pow_even`: $(-x)^{2k} = x^{2k}$ on `Rat`.
* `neg_pow_odd`: $(-x)^{2k+1} = - x^{2k+1}$ on `Rat`.
* `even_pole_pair_doubles`: $\frac{1}{n^{2k}} + \frac{1}{(-n)^{2k}} = \frac{2}{n^{2k}}$
  (why Euler's cotangent method extracts $2\zeta(2k) \in \mathbb{Q}\cdot \pi^{2k}$ via Bernoulli numbers $B_{2k}$).
* `odd_pole_pair_cancels`: $\frac{1}{n^{2k+1}} + \frac{1}{(-n)^{2k+1}} = 0$
  (why the exact same method yields $0 = 0$ for $\zeta(3), \zeta(5), \dots$).
* `symmetric_partial_sum_even` and `symmetric_partial_sum_odd`: the finite symmetric partial-sum versions.
-/

namespace Zeta3

theorem neg_pow_even (x : Rat) : ∀ k : Nat, (-x) ^ (2 * k) = x ^ (2 * k) := by
  intro k
  induction k with
  | zero =>
      simp
  | succ m ih =>
      have h1 : 2 * (m + 1) = 2 * m + 2 := by omega
      rw [h1]
      have hpow : ∀ (a : Rat) (n : Nat), a ^ (n + 2) = a ^ n * (a * a) := by
        intro a n
        grind
      rw [hpow (-x) (2 * m), hpow x (2 * m), ih]
      grind

theorem neg_pow_odd (x : Rat) (k : Nat) : (-x) ^ (2 * k + 1) = -(x ^ (2 * k + 1)) := by
  have hpow : ∀ (a : Rat) (n : Nat), a ^ (n + 1) = a ^ n * a := by
    intro a n
    grind
  rw [hpow (-x) (2 * k), hpow x (2 * k), neg_pow_even x k]
  grind

/-- Symmetric pole-pairing at `+n` and `-n` **doubles** for even powers `2k` (e.g. `ζ(2), ζ(4), …`). -/
theorem even_pole_pair_doubles (x : Rat) (k : Nat) :
    (1 : Rat) / (x ^ (2 * k)) + (1 : Rat) / ((-x) ^ (2 * k)) = (2 : Rat) / (x ^ (2 * k)) := by
  rw [neg_pow_even x k]
  grind

/-- Helper: `(-a)⁻¹ = -(a⁻¹)` on `Rat`. -/
theorem rat_inv_neg (a : Rat) : (-a)⁻¹ = -(a⁻¹) := by
  by_cases ha : a = 0
  · subst ha; decide +kernel
  · have hna : -a ≠ 0 := by grind
    have hmul : (-a) * (-(a⁻¹)) = 1 := by
      have h1 : (-a) * (-(a⁻¹)) = a * a⁻¹ := by grind
      rw [h1, Rat.mul_inv_cancel a ha]
    have hmul2 : (-a) * (-a)⁻¹ = 1 := Rat.mul_inv_cancel (-a) hna
    have heq : (-a) * (-a)⁻¹ = (-a) * (-(a⁻¹)) := by rw [hmul2, hmul]
    have hleft := congrArg (fun z => (-a)⁻¹ * z) heq
    rw [← Rat.mul_assoc, ← Rat.mul_assoc, Rat.inv_mul_cancel (-a) hna,
        Rat.one_mul, Rat.one_mul] at hleft
    exact hleft

/-- Symmetric pole-pairing at `+n` and `-n` **cancels identically** for odd powers `2k + 1` (e.g. `ζ(3), ζ(5), …`). -/
theorem odd_pole_pair_cancels (x : Rat) (k : Nat) :
    (1 : Rat) / (x ^ (2 * k + 1)) + (1 : Rat) / ((-x) ^ (2 * k + 1)) = 0 := by
  rw [neg_pow_odd x k]
  have hdiv : ∀ (a : Rat), (1 : Rat) / (-a) = -((1 : Rat) / a) := by
    intro a
    rw [Rat.div_def, Rat.div_def, Rat.one_mul, Rat.one_mul]
    exact rat_inv_neg a
  rw [hdiv]
  grind

/-- Helper: partial sum of a function `f : Nat → Rat` over `n = 1, …, N`. -/
def partialSum (f : Nat → Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => partialSum f n + f (n + 1)

/-- For every truncation `N` and every even exponent `2k`, the symmetric sum doubles `∑_{n=1}^N 1/n^{2k}`. -/
theorem symmetric_partial_sum_even (k : Nat) :
    ∀ N : Nat,
      partialSum (fun n => (1 : Rat) / ((n : Rat) ^ (2 * k)) + (1 : Rat) / ((-(n : Rat)) ^ (2 * k))) N =
      2 * partialSum (fun n => (1 : Rat) / ((n : Rat) ^ (2 * k))) N := by
  intro N
  induction N with
  | zero =>
      simp [partialSum]
  | succ m ih =>
      simp only [partialSum, ih]
      rw [even_pole_pair_doubles ((m + 1 : Nat) : Rat) k]
      grind

/-- For every truncation `N` and every odd exponent `2k+1` (including `ζ(3)` at `k=1`), the symmetric sum vanishes identically. -/
theorem symmetric_partial_sum_odd (k : Nat) :
    ∀ N : Nat,
      partialSum (fun n => (1 : Rat) / ((n : Rat) ^ (2 * k + 1)) + (1 : Rat) / ((-(n : Rat)) ^ (2 * k + 1))) N = 0 := by
  intro N
  induction N with
  | zero =>
      simp [partialSum]
  | succ m ih =>
      simp only [partialSum, ih]
      rw [odd_pole_pair_cancels ((m + 1 : Nat) : Rat) k]
      grind

/-- Specialized to `ζ(3)` (`k = 1`): symmetric residue pairing vanishes identically for every truncation `N`. -/
theorem zeta3_symmetric_residue_cancels (N : Nat) :
    partialSum (fun n => (1 : Rat) / ((n : Rat) ^ 3) + (1 : Rat) / ((-(n : Rat)) ^ 3)) N = 0 :=
  symmetric_partial_sum_odd 1 N

end Zeta3
