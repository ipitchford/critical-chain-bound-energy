"""Rank-match critical Ising annulus spectra with exact integer multiplicities.

Free/free gives the Neveu–Schwarz product over positive half-integers, and
free/fixed gives the Ramond product over positive integers. All reported
energies are measured above their individual ground levels in units pi*v/L.
"""
import csv
import math
import sys


def require(condition, message):
    """Explicit check that survives python -O (no assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED: " + message)



def euler_terms(scale, cutoff):
    """Terms of Euler's pentagonal series other than the constant."""
    k = 1
    while True:
        small = scale * k * (3 * k - 1) // 2
        if small > cutoff:
            break
        sign = -1 if k % 2 else 1
        yield small, sign
        large = scale * k * (3 * k + 1) // 2
        if large <= cutoff:
            yield large, sign
        k += 1


def euler_product(series, scale, invert=False):
    """Multiply by prod(1-q^(scale*j)) or divide as a formal series."""
    result = series.copy()
    cutoff = len(series) - 1
    if invert:
        terms = list(euler_terms(scale, cutoff))
        for k in range(1, cutoff + 1):
            value = result[k]
            for offset, sign in terms:
                if offset > k:
                    break
                value -= sign * result[k - offset]
            result[k] = value
    else:
        for offset, sign in euler_terms(scale, cutoff):
            for k in range(offset, cutoff + 1):
                result[k] += sign * series[k - offset]
    return result


def multiplicities(cutoff):
    # On x=q^(1/2) grid, NS product is E(x^2)^2/(E(x)E(x^4)),
    # while R counted on *integer* grid is E(q^2)/E(q).
    seed = [1] + [0] * (2 * cutoff)
    seed = euler_product(seed, 2)
    ramond = euler_product(seed, 1, invert=True)
    ns = euler_product(ramond, 2)
    ns = euler_product(ns, 4, invert=True)
    require(ns[:10] == [1, 1, 0, 1, 1, 1, 1, 1, 2, 2], "check_boundary_characters.py:55")
    require(ramond[:10] == [1, 1, 1, 2, 2, 3, 4, 5, 6, 8], "check_boundary_characters.py:56")
    return ns, ramond


def matched_energies(ns, ramond):
    level_r = 0
    available = ramond[0]
    summed_energies = []
    for multiplicity in ns:
        needed = multiplicity
        energy_sum = 0
        while needed:
            if not available:
                level_r += 1
                available = ramond[level_r]
                continue
            block = min(needed, available)
            energy_sum += block * level_r
            needed -= block
            available -= block
        summed_energies.append(energy_sum)
    return summed_energies


def free_stats(u):
    # With u≥0.008 in the supplied runs, 140 decay lengths give negligible error.
    cutoff = math.ceil(140 / u)
    levels = [j + 0.5 for j in range(cutoff)]
    log_z = math.fsum(math.log1p(math.exp(-u * h)) for h in levels)
    mean = math.fsum(h / (1 + math.exp(u * h)) for h in levels)
    return log_z, mean


def calculate():
    cutoff = 13000
    ns, ramond = multiplicities(cutoff)
    matched = matched_energies(ns, ramond)
    K = math.pi**2 / 12
    boundary_shift = math.log(2) / 2
    for u in (0.20, 0.10, 0.06, 0.04, 0.025, 0.016):
        log_z, mean = free_stats(u)
        S = log_z + u * mean
        beta = math.sqrt(K / cutoff)
        log_z_tilted, mean_tilted = free_stats(beta)
        # Map occupied NS modes j+1/2 to R modes j+1. This is a bijection
        # with configurationwise R energy <= twice NS energy, so this bound
        # also holds pointwise for sorted spectra and bounds the omitted tail.
        log_tail = math.log(2 * mean_tilted) - (u - beta) * cutoff + log_z_tilted - log_z
        passive = math.fsum(x * math.exp(-u * j / 2 - log_z) for j, x in enumerate(matched))
        mass = math.fsum(n * math.exp(-u * j / 2 - log_z) for j, n in enumerate(ns))
        require(abs(mass - 1) < 1e-9, "check_boundary_characters.py:106")
        yield dict(u=u, entropy=S, passive=passive, recovered_mass=mass,
                   log10_truncation_energy_bound=log_tail / math.log(10),
                   inferred_boundary_shift=2 * math.sqrt(K * passive) - S,
                   predicted_boundary_shift=boundary_shift)


if __name__ == '__main__':
    rows = list(calculate())
    writer = csv.DictWriter(sys.stdout, rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
