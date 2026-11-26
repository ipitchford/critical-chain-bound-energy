"""Finite-system checks for passive energies of odd-power open chains.

Run: python check_passive_fermi.py > passive_fermi.csv
Requires Python >=3.10, NumPy and SciPy. Values are in unit hopping energy.
The unrestricted minimum orders all Fock-space eigenvalues, including every
particle-number sector. The alternative Gaussian minimum orders one-body modes.
"""
import csv
import itertools
import sys

import numpy as np
from scipy.special import entr


def require(condition, message):
    """Explicit check that survives python -O (no assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED: " + message)



def path_adjacency(n):
    return np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)


def subset_sums(cost):
    totals = np.array([0.0])
    for value in cost:
        totals = np.concatenate([totals, totals + value])
    return np.sort(totals)


def fock_probabilities(occupations):
    probabilities = np.array([1.0])
    for nu in occupations:
        probabilities = np.concatenate(
            [probabilities * (1 - nu), probabilities * nu]
        )
    return np.sort(probabilities)[::-1]


def sector_vectors(values, as_occupations):
    sectors = [np.array([1.0 if as_occupations else 0.0])]
    for value in values:
        if as_occupations:
            empty = [part * (1 - value) for part in sectors]
            full = [part * value for part in sectors]
        else:
            empty = sectors
            full = [part + value for part in sectors]
        sectors = [empty[0]] + [
            np.concatenate([empty[n], full[n - 1]]) for n in range(1, len(sectors))
        ] + [full[-1]]
    return [np.sort(part) for part in sectors]


def exact_small_rdm(length, occupied):
    orbitals = occupied.shape[0]
    vector = np.zeros(1 << orbitals)
    for sites in itertools.combinations(range(orbitals), orbitals // 2):
        bits = sum(1 << site for site in sites)
        vector[bits] = np.linalg.det(occupied[np.array(sites), :])
    require(abs(np.dot(vector, vector) - 1) < 1e-11, "check_passive_fermi.py:57")
    matrix = vector.reshape(1 << length, 1 << length).T
    return matrix @ matrix.T


def many_body_hamiltonian(h):
    """Direct Jordan–Wigner matrix using bit operations, independent of subset sums."""
    modes = len(h)
    ham = np.zeros((1 << modes, 1 << modes))
    for row_site, column_site in np.argwhere(h != 0):
        i, j = int(row_site), int(column_site)
        for ket in range(1 << modes):
            if not (ket & (1 << j)):
                continue
            tmp = ket ^ (1 << j)
            if tmp & (1 << i):
                continue
            parity = (ket & ((1 << j) - 1)).bit_count()
            parity += (tmp & ((1 << i) - 1)).bit_count()
            bra = tmp | (1 << i)
            ham[bra, ket] += (-1) ** parity * h[i, j]
    return ham


def check(length, writer):
    total = 2 * length
    path = path_adjacency(total)
    energies, modes = np.linalg.eigh(path)
    occupied = modes[:, energies < 0]
    correlation = occupied[:length] @ occupied[:length].T
    nu = np.clip(np.linalg.eigvalsh(correlation), 0.0, 1.0)
    entropy = np.sum(entr(nu) + entr(1 - nu))
    weights = fock_probabilities(nu)
    sector_weights = sector_vectors(nu, as_occupations=True)
    require(abs(np.sum(weights) - 1) < 1e-12, "check_passive_fermi.py:91")
    if length <= 6:
        state = exact_small_rdm(length, occupied)
        state_weights = np.linalg.eigvalsh(state)[::-1]
        require(np.max(np.abs(state_weights - weights)) < 1e-12, "check_passive_fermi.py:95")
    for z in (1, 3, 5):
        block = np.linalg.matrix_power(path, z)[:length, :length]
        epsilon = np.linalg.eigvalsh(block)
        excitation = subset_sums(abs(epsilon))
        passive = np.dot(weights, excitation)
        sector_energies = sector_vectors(epsilon, as_occupations=False)
        gs = np.sum(epsilon[epsilon < 0])
        conserving = sum(
            np.dot(p[::-1], e - gs) for p, e in zip(sector_weights, sector_energies)
        )
        # For each aligned signed one-body energy, deviations from the
        # ground occupation are nu for positive levels and 1-nu otherwise.
        aligned = nu[::-1]
        mismatch = np.where(epsilon >= 0, aligned, 1 - aligned)
        gaussian = np.dot(abs(epsilon), mismatch)
        require(gaussian >= conserving - 2e-12, "check_passive_fermi.py:111")
        require(conserving >= passive - 2e-12, "check_passive_fermi.py:112")
        direct = ""
        if length <= 6:
            many_energies = np.linalg.eigvalsh(many_body_hamiltonian(block))
            many_energies -= many_energies[0]
            require(np.max(abs(many_energies - excitation)) < 2e-12, "check_passive_fermi.py:117")
            direct = np.dot(state_weights, many_energies)
            require(abs(direct - passive) < 2e-12, "check_passive_fermi.py:119")
        writer.writerow(
            [length, total, z, f"{entropy:.13g}", f"{passive:.13g}",
             f"{conserving:.13g}", f"{gaussian:.13g}", f"{gaussian-passive:.9g}",
             "" if direct == "" else f"{direct:.13g}"]
        )


if __name__ == "__main__":
    output = csv.writer(sys.stdout, lineterminator="\n")
    output.writerow(["L", "N", "z", "entropy_nats", "Q_full", "Q_U1", "Q_one_body",
                     "one_body_minus_full", "Q_direct_fock"])
    for length in (2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16):
        check(length, output)
