#!/usr/bin/env python3
"""Rank-truncated passive-energy benchmark for the open critical Ising chain.

Python 3, NumPy and SciPy. For unrounded input, retained ranks give a lower
bound on the passive excess and the unretained mass times spectral width
bounds the truncation error. These are floating-point estimates, not intervals
certified against rounding. Compare the first entries with reproduce_ising.py.
"""
import csv
import math
from pathlib import Path

import numpy as np
from scipy.linalg import svdvals


def restricted_cross(size, r):
    """Sine-projector Toeplitz kernel at odd integer arguments r."""
    signs = np.where(r % 4 == 1, 1, -1)
    return signs / (2 * (size + 1) * np.sin(np.pi * r / (2 * (size + 1))))


def ising_excitation_probabilities(length):
    # Γ_AB is 2 C_AB up to phases. One parity block contains one
    # singular value from each pair associated to one entangled Ising mode.
    size = 4 * length
    inside = np.arange(1, 2 * length + 1, 2)[:, None]
    outside = np.arange(2 * length + 2, size + 1, 2)[None, :]
    cross = restricted_cross(size, inside - outside) - restricted_cross(size, inside + outside)
    singular = svdvals(cross, overwrite_a=True)
    t = singular**2
    # Stable root of p(1-p)=t. All p are the smaller mode occupations.
    return 2 * t / (1 + np.sqrt(np.maximum(0.0, 1 - 4 * t)))


def lowest_subset_sums(costs, retain):
    """Keep the first retain subset sums with multiplicities using positivity."""
    sums = np.array([0.0])
    for energy in np.sort(costs):
        if len(sums) == retain and energy >= sums[-1]:
            break
        expanded = np.concatenate((sums, energy + sums))
        if len(expanded) > retain:
            expanded = np.partition(expanded, retain - 1)[:retain]
        sums = np.sort(expanded)
    return sums


def estimate(length, retain):
    p = ising_excitation_probabilities(length)
    p = p[p > 0]
    entropic = np.log1p(-p) - np.log(p)
    entropy = math.fsum(-p*np.log(p) - (1-p)*np.log1p(-p))
    log_pmax = math.fsum(np.log1p(-p))
    weights = np.exp(log_pmax - lowest_subset_sums(entropic, retain))
    j = np.arange(length, dtype=float)
    costs = 4 * np.sin(np.pi*(2*j+1)/(4*length+2))
    energies = lowest_subset_sums(costs, retain)
    missing = max(0.0, 1.0-math.fsum(weights))
    # Absolute, deliberately conservative allowance for the scalar mass
    # subtraction. This does not certify the matrix SVD against rounding.
    tail_estimate = (missing + 8*np.finfo(float).eps) * math.fsum(costs)
    partial = weights @ energies
    denom = entropy**2/length
    return [length, len(weights), entropy, partial, missing, tail_estimate,
            partial/denom, tail_estimate/denom]


if __name__ == '__main__':
    fields = ['block_length', 'retained_ranks', 'entropy', 'partial_passive_excess',
              'estimated_missing_mass', 'tail_energy_envelope', 'L_Qpartial_over_S2',
              'tail_envelope_times_L_over_S2']
    with (Path(__file__).resolve().parent.parent / 'results' / 'ising_mesoscopic.csv').open('w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(fields)
        for n in (4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048):
            row = estimate(n, 16384)
            writer.writerow(row)
            file.flush()
            print(f'L={n:5d}, L Q / S²={row[-2]:.8f}, tail allowance={row[-1]:.2e}', flush=True)
