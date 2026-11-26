#!/usr/bin/env python3
"""Check the critical open Ising passive energy in spin and Majorana representations.

Run with Python 3, NumPy and SciPy. The output file is an optional argument.
The default is ising_checks.csv next to this script. L is the block length.
The occupation-number basis contains every local parity and charge sector.
The finite-size ratios test normalisation; they are not asymptotic fits.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh


def require(condition, message):
    """Explicit check that survives python -O (no assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED: " + message)



def spin_chain(n):
    """Open Pauli Hamiltonian -sum(XX)-sum(Z), in the computational basis."""
    dim = 2**n
    indices = np.arange(dim)
    diagonal = np.fromiter((2 * index.bit_count() - n for index in range(dim)),
                           dtype=float, count=dim)
    rows = [indices]
    columns = [indices]
    entries = [diagonal]
    for site in range(n - 1):
        rows.append(indices)
        columns.append(indices ^ (3 << site))
        entries.append(-np.ones(dim))
    return sparse.coo_matrix((np.concatenate(entries),
                              (np.concatenate(rows), np.concatenate(columns))),
                             shape=(dim, dim)).tocsr()


def subset_sums(costs):
    sums = np.array([0.0])
    for cost in costs:
        sums = np.r_[sums, sums + cost]
    return np.sort(sums)


def binary_spectrum(probabilities):
    weights = np.array([1.0])
    for p in probabilities:
        weights = np.r_[weights * p, weights * (1 - p)]
    return np.sort(weights)[::-1]


def entropy(weights):
    positive = weights[weights > 0]
    return -positive @ np.log(positive)


def check_block(length):
    dim = 2**length
    spin = spin_chain(2 * length)
    generator = np.random.default_rng(23023 + length)
    _, vector = eigsh(spin, k=1, which='SA', tol=2e-14,
                      v0=generator.normal(size=dim**2))
    singular = np.linalg.svd(vector[:, 0].reshape(dim, dim), compute_uv=False)
    spin_weights = singular**2
    block_energies = np.linalg.eigvalsh(spin_chain(length).toarray())
    spin_costs = block_energies - block_energies[0]

    size = 4 * length
    occupied_indices = np.arange(size // 2 + 1, size + 1)
    block_indices = np.arange(1, size // 2 + 1)
    occupied = np.sqrt(2.0 / (size + 1)) * np.sin(
        np.pi * np.outer(block_indices, occupied_indices) / (size + 1))
    xx_occupation = np.clip(np.linalg.eigvalsh(occupied @ occupied.T), 0, 1)
    # Each Majorana pair gives one Ising occupation and two XX occupations.
    ising_weights = binary_spectrum(xx_occupation[length:])
    phase = np.arange(1, length + 1) * np.pi / (2 * length + 1)
    ising_costs = subset_sums(4 * np.cos(phase))

    spin_q = spin_weights @ spin_costs
    ising_q = ising_weights @ ising_costs
    spin_entropy = entropy(spin_weights)
    xx_entropy = np.sum([entropy(np.array([p, 1 - p])) for p in xx_occupation])
    spectrum_error = max(abs(spin_weights - ising_weights))
    energy_error = max(abs(spin_costs - ising_costs))
    require(spectrum_error < 2e-13 and energy_error < 2e-12, "reproduce_ising.py:83")
    require(abs(spin_entropy - xx_entropy / 2) < 2e-12, "reproduce_ising.py:84")
    require(abs(spin_q - ising_q) < 2e-13, "reproduce_ising.py:85")
    renyi_errors = []
    for order in (0.5, 1.5, 2., 3.):
        ising_logtrace = np.log(np.sum(spin_weights**order))
        xx_logtrace = sum(np.log(p**order + (1-p)**order) for p in xx_occupation)
        renyi_errors.append(abs(ising_logtrace - xx_logtrace / 2))
    # At q=1/2, square roots magnify roundoff at boundary occupations.
    require(max(renyi_errors) < 2e-7, "reproduce_ising.py:92")
    return dict(block_length=length, full_spins=2*length,
                entropy_spin=spin_entropy, half_entropy_xx=xx_entropy/2,
                passive_excess_spin=spin_q, passive_excess_majorana=ising_q,
                block_times_excess_over_entropy_squared=length*spin_q/spin_entropy**2,
                max_schmidt_probability_error=spectrum_error,
                max_block_energy_error=energy_error,
                max_logtrace_error= max(renyi_errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', nargs='?', type=Path,
                        default=Path(__file__).resolve().parent.parent / 'results' / 'ising_checks.csv')
    parser.add_argument('--max-length', type=int, default=8)
    args = parser.parse_args()
    results = [check_block(n) for n in range(1, args.max_length + 1)]
    with args.output.open('w', newline='') as dest:
        writer = csv.DictWriter(dest, fieldnames=results[0])
        writer.writeheader()
        writer.writerows(results)
    print(f"{len(results)} sizes verified. Final ratio {results[-1]['block_times_excess_over_entropy_squared']:.8g}.")
    print(f"Asymptotic prediction for this normalisation: {6 / np.pi:.8g}")


if __name__ == '__main__':
    main()
