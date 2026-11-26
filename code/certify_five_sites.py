"""Interval check of the ordering inversion at N=10, L=5, z=1.

Uses mpmath's interval arithmetic. Every comparison uses separated intervals.
The additive energy unit is the path hopping. No numerical eigensolver is used.
"""
from mpmath import iv


def require(condition, message):
    """Explicit check that survives python -O (no assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED: " + message)


iv.dps = 40
half = iv.mpf(1) / 2


def kernel(r):
    if r == 0:
        return half
    # r is an integer: evaluate sin(pi*r/2) as its exact value.
    sine_cycle = (0, 1, 0, -1)
    return sine_cycle[r % 4] / (22 * iv.sin(iv.pi * r / 22))


# In odd/even site order the off-diagonal block is 3 x 2.
cross = [[kernel(i - j) - kernel(i + j) for j in (2, 4)] for i in (1, 3, 5)]
a = sum(row[0] ** 2 for row in cross)
d = sum(row[1] ** 2 for row in cross)
b = sum(row[0] * row[1] for row in cross)
discriminant = iv.sqrt((a - d) ** 2 + 4 * b ** 2)
nu_far = half - iv.sqrt((a + d + discriminant) / 2)
nu_near = half - iv.sqrt((a + d - discriminant) / 2)
require(nu_far > 0 and nu_near > nu_far and nu_near < half, "certify_five_sites.py:28")
far_odds = nu_far / (1 - nu_far)
near_odds = nu_near / (1 - nu_near)
probability_ratio = near_odds ** 2 / far_odds
energy_inversion = 2 - iv.sqrt(3)
require(energy_inversion > 0 and probability_ratio > 1, "certify_five_sites.py:33")
print("weight(two near flips)/weight(one far flip) =", probability_ratio)
print("energy(two near flips)-energy(one far flip) =", energy_inversion)
