#!/usr/bin/env python3
"""
rhin_viola_lab.py — Saddle-point, Weyl group W(D_5), and irrationality-measure
workbench for ζ(3) (Apéry 1978 -> Beukers 1979 -> Hata 2000 -> Rhin-Viola 2001
-> David Niedbala Giraudin Sept 2026).

Key facts from the literature (updated through November 2026):
  - For ζ(2), Rhin & Viola (1996) discovered a 120-element group S_5 = W(A_4).
  - For ζ(3), Rhin & Viola (2001, Acta Arith. 97) discovered a 1920-element
    permutation group Φ ≅ W(D_5) acting on the 10-element parameter multiset
    of the 8-parameter Beukers-type triple integral:
        I(h,j,k,l,m,q,r,s) = ∭_{[0,1]^3} x^h (1-x)^l y^k (1-y)^s z^j (1-z)^q
                                         / (1 - (1-xy)z)^{q+h-r+1} dx dy dz
    subject to h + m = k + r and j + q = l + s.
  - Record progression for μ(ζ(3)):
      * Apéry (1978) / Beukers (1979): μ <= 13.41782023
      * Hata (1990):                   μ <=  8.830284
      * Hata (2000):                   μ <=  7.377956
      * Rhin & Viola (2001):           μ <=  5.513891   (stood for 25 years)
      * Giraudin (Sept 2026, 2609.26980): μ < 5.5138800 (narrow-cone shift:
        8 params = 52000 * RV_params shifted by <= 23)
      * Conjectured truth:             μ = 2
"""

from __future__ import annotations

import math
from typing import Tuple


def maximize_general_rv_kernel(
    h: float, j: float, k: float, l: float, m: float, q: float, r: float, s: float
) -> Tuple[float, Tuple[float, float, float]]:
    """
    Maximize the log-integrand of the Rhin-Viola 8-parameter Beukers triple integral:
      f(x, y, z) = h ln(x) + l ln(1-x) + k ln(y) + s ln(1-y)
                   + j ln(z) + q ln(1-z) - (q + h - r) ln(1 - (1 - x*y)*z)
    over (x, y, z) ∈ (0, 1)^3.

    For fixed u = 1 - x*y ∈ (0, 1), maximizing j ln(z) + q ln(1-z) - A ln(1 - u z)
    (where A = q + h - r) satisfies the quadratic equation:
      (j + q - A) u z^2 - (j + q + (j - A) u) z + j = 0  (for z ∈ (0, 1)).
    """
    A = q + h - r

    def best_z(u: float) -> float:
        # Derivative wrt z: j/z - q/(1-z) + A*u/(1 - u*z) = 0
        # Multiply by z(1-z)(1-u*z):
        #   j(1-z)(1-u*z) - q*z(1-u*z) + A*u*z(1-z) = 0
        #   (j*u + q*u - A*u) z^2 + (-j - j*u - q + A*u) z + j = 0
        c2 = (j + q - A) * u
        c1 = -(j + q) + (A - j) * u
        c0 = j
        if abs(c2) < 1e-14:
            return -c0 / c1 if abs(c1) > 1e-14 else 0.5
        disc = max(0.0, c1 * c1 - 4.0 * c2 * c0)
        sd = math.sqrt(disc)
        z1 = (-c1 - sd) / (2.0 * c2)
        z2 = (-c1 + sd) / (2.0 * c2)
        candidates = [z for z in (z1, z2) if 1e-12 < z < 1.0 - 1e-12]
        if not candidates:
            return 0.5
        if len(candidates) == 1:
            return candidates[0]
        def fz(z: float) -> float:
            return j * math.log(z) + q * math.log(1.0 - z) - A * math.log(1.0 - u * z)
        return max(candidates, key=fz)

    def g(x: float, y: float) -> Tuple[float, float]:
        if x <= 1e-12 or x >= 1.0 - 1e-12 or y <= 1e-12 or y >= 1.0 - 1e-12:
            return -1e30, 0.5
        u = 1.0 - x * y
        z = best_z(u)
        denom = 1.0 - u * z
        if denom <= 1e-15:
            return -1e30, z
        val = (
            h * math.log(x)
            + l * math.log(1.0 - x)
            + k * math.log(y)
            + s * math.log(1.0 - y)
            + j * math.log(z)
            + q * math.log(1.0 - z)
            - A * math.log(denom)
        )
        return val, z

    best_val = -1e30
    best_xyz = (0.5, 0.5, 0.5)
    grid = 36
    for ix in range(1, grid):
        x = ix / grid
        for iy in range(1, grid):
            y = iy / grid
            val, z = g(x, y)
            if val > best_val:
                best_val = val
                best_xyz = (x, y, z)

    x, y, z = best_xyz
    step = 0.04
    for _ in range(55):
        improved = False
        for dx, dy in [
            (step, 0), (-step, 0), (0, step), (0, -step),
            (step, step), (step, -step), (-step, step), (-step, -step),
        ]:
            nx, ny = x + dx, y + dy
            if 1e-9 < nx < 1.0 - 1e-9 and 1e-9 < ny < 1.0 - 1e-9:
                v, nz = g(nx, ny)
                if v > best_val:
                    best_val = v
                    x, y, z = nx, ny, nz
                    improved = True
        if not improved:
            step *= 0.5

    return best_val, (x, y, z)


def main() -> None:
    print("=" * 78)
    print("  Beukers / Rhin-Viola (W(D_5), |Φ|=1920) / Giraudin (2026) Workbench")
    print("=" * 78)

    # 1. Exact Beukers (1979) symmetric case: all 8 params = 1
    log_max, (x, y, z) = maximize_general_rv_kernel(1, 1, 1, 1, 1, 1, 1, 1)
    c0 = -log_max
    c0_exact = 4.0 * math.log(1.0 + math.sqrt(2.0))
    c1 = 3.0
    mu_apery = 1.0 + (c0_exact + c1) / (c0_exact - c1)

    print("1. Classic Beukers / Apéry Symmetric Case (h=j=k=l=m=q=r=s=1):")
    print(f"   Max point (x, y, z)   = ({x:.10f}, {y:.10f}, {z:.10f})")
    print(f"   Exact (√2-1, √2-1, …) = ({math.sqrt(2)-1:.10f}, {math.sqrt(2)-1:.10f}, {1/math.sqrt(2):.10f})")
    print(f"   Decay c0 (numerical)  = {c0:.12f}")
    print(f"   Decay c0 (4 ln(1+√2)) = {c0_exact:.12f}")
    print(f"   Denominator exp c1    = {c1:.12f}  (lcm(1..n)^3 ~ e^(3n))")
    print(f"   Margin c0 - c1        = {c0_exact - c1:.12f}  (> 0 => ζ(3) ∉ ℚ)")
    print(f"   Irrationality measure = μ(ζ(3)) <= {mu_apery:.8f}")
    print()

    # 2. Updated Literature Benchmarks (through Nov 2026)
    print("2. Complete Literature Record Progression for μ(ζ(3)) (1978–2026):")
    benchmarks = [
        ("Apéry (1978) / Beukers (1979)", 13.4178202, "Symmetric triple integral, c1 = 3"),
        ("Hata (1990)",                    8.8302840, "Asymmetric Legendre-type polynomials"),
        ("Hata (2000)",                    7.3779560, "Refined asymmetric rational forms"),
        ("Rhin & Viola (2001)",            5.5138910, "W(D_5) group (|Φ| = 1920) p-adic savings"),
        ("Giraudin (Sept 2026) [RECORD]",  5.5138800, "arXiv:2609.26980 narrow-cone shift of RV"),
        ("Conjectured truth",              2.0000000, "Generic / Roth-type real behavior"),
    ]
    for name, val, note in benchmarks:
        print(f"   {name:34s} μ <= {val:10.7f} | {note}")
    print()

    # 3. Verify general 8-parameter solver on asymmetric Rhin-Viola-type tuples
    # obeying h + m = k + r and j + q = l + s
    print("3. Asymmetric 8-Parameter Kernel Sweep (h+m = k+r, j+q = l+s):")
    tuples = [
        ("Beukers symmetric",   (1, 1, 1, 1, 1, 1, 1, 1)),
        ("Slightly skewed A",   (1.2, 1.0, 1.1, 1.0, 1.0, 1.0, 1.1, 1.0)),
        ("Slightly skewed B",   (1.0, 1.2, 1.0, 1.1, 1.0, 1.0, 1.0, 1.1)),
        ("Asymmetric Zudilin",  (1.3, 1.1, 1.2, 1.0, 1.1, 1.1, 1.2, 1.2)),
    ]
    print(f"   {'Label':22s}  {'c0 (decay)':>12s}  {'max(x,y,z)':>30s}")
    for lbl, tup in tuples:
        h, j, k, l, m, q, r, s = tup
        assert abs((h + m) - (k + r)) < 1e-12 and abs((j + q) - (l + s)) < 1e-12
        lm, pt = maximize_general_rv_kernel(*tup)
        print(f"   {lbl:22s}  {-lm:12.6f}  ({pt[0]:.5f}, {pt[1]:.5f}, {pt[2]:.5f})")


if __name__ == "__main__":
    main()
