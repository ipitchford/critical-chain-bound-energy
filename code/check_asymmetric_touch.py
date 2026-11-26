"""Finite-system checks for the range-three one-sided touch, a = 3/2.

Run python code/check_asymmetric_touch.py > results/asymmetric_touch.csv; it imports
the sibling module check_passive_fermi.py. All unitaries called G-number conserve particle number.
For L <= 6, explicit Slater partial trace and direct Fock matrices independently
verify both unrestricted and number-conserving residuals.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # sibling module, works under python -I

import numpy as np
from scipy.special import entr

from check_passive_fermi import (
    exact_small_rdm,
    fock_probabilities,
    many_body_hamiltonian,
    path_adjacency,
    sector_vectors,
    subset_sums,
)


def require(condition, message):
    """Explicit check that survives python -O (no assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED: " + message)



def measure(length):
    total = 2 * length
    path = path_adjacency(total)
    numbers, modes = np.linalg.eigh(path)
    occupied = modes[:, numbers < 0]
    nu = np.clip(np.linalg.eigvalsh(occupied[:length] @ occupied[:length].T), 0, 1)
    fermi_weight = fock_probabilities(nu)
    a = 1.5
    shifted = path - a * np.eye(total)
    full = path @ shifted @ shifted
    require(np.linalg.norm(full @ occupied - occupied * (numbers[numbers < 0] * (numbers[numbers < 0] - a) ** 2)) < 2e-12, "check_asymmetric_touch.py:35")
    block = full[:length, :length]
    eps = np.linalg.eigvalsh(block)
    gs = np.sum(eps[eps < 0])
    q_all = np.dot(fermi_weight, subset_sums(abs(eps)))
    weights = sector_vectors(nu, as_occupations=True)
    energies = sector_vectors(eps, as_occupations=False)
    q_number = sum(np.dot(w[::-1], e - gs) for w, e in zip(weights, energies))
    # Finite charge-floor bound uses only the weight of each sector.
    q_floor = sum(np.sum(w) * (e[0] - gs) for w, e in zip(weights, energies))
    penetration = np.where(eps >= 0, nu[::-1], 1 - nu[::-1])
    q_gauss = np.dot(abs(eps), penetration)
    require(q_number >= max(q_all, q_floor) - 2e-12, "check_asymmetric_touch.py:47")
    require(q_gauss >= q_number - 2e-12, "check_asymmetric_touch.py:48")
    independent_all = ""
    independent_number = ""
    if length <= 6:
        state = exact_small_rdm(length, occupied)
        H = many_body_hamiltonian(block)
        energies_all = np.linalg.eigvalsh(H)
        minimum = energies_all[0]
        probabilities = np.linalg.eigvalsh(state)[::-1]
        direct_all = np.dot(probabilities, energies_all - minimum)
        direct_number = 0.0
        for n in range(length + 1):
            selected = [m for m in range(1 << length) if m.bit_count() == n]
            ix = np.ix_(selected, selected)
            probabilities = np.linalg.eigvalsh(state[ix])[::-1]
            energies_n = np.linalg.eigvalsh(H[ix]) - minimum
            direct_number += np.dot(probabilities, energies_n)
        require(abs(direct_all - q_all) < 2e-12, "check_asymmetric_touch.py:65")
        require(abs(direct_number - q_number) < 2e-12, "check_asymmetric_touch.py:66")
        independent_all = f"{direct_all:.13g}"
        independent_number = f"{direct_number:.13g}"
    return [length, a, f"{np.sum(entr(nu) + entr(1-nu)):.13g}",
            f"{np.sum(nu * (1-nu)):.13g}", np.sum(eps < -1e-12),
            *(f"{v:.13g}" for v in (q_all, q_floor, q_number, q_gauss)),
            independent_all, independent_number]


if __name__ == "__main__":
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["L", "a", "entropy_nats", "VarN", "n_negative",
                     "Q_full", "Q_charge_floor", "Q_U1", "Q_G_number",
                     "direct_Q_full", "direct_Q_U1"])
    for length in (2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16):
        writer.writerow(measure(length))
