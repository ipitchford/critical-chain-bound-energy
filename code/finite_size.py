"""Finite-size values of the residual energies against their L -> infinity limits.

For an open hopping chain of 2L sites, the half-filled ground state of dGamma(f(T_2L)) is the
same for every sign-preserving polynomial f.  This script computes, for the left block of L sites:

  S      block entropy (nats)
  Q      unrestricted residual: sum_r p_r E_r with both lists ranked (truncated ranking, with a
         certified bound on the omitted probability mass and its energy contribution)
  QG     number-conserving Gaussian residual (signed one-body rearrangement)
  floor  a rigorous lower bound on the charge-conserving residual Q^ch: the expected cost of the
         holes (and particles) that the block charge n forces, sum of the cheapest |n_* - n| levels

and the limits proved in the manuscript.  Every check raises SystemExit on failure (no assert).
Usage: python3 code/finite_size.py > results/finite_size.csv
"""
import argparse
import heapq
import math
import sys

import mpmath as mp
import numpy as np
from scipy import sparse
from scipy.linalg import eigvalsh

mp.mp.dps = 30
LOG2 = math.log(2.0)


def die(msg):
    sys.stderr.write("CHECK FAILED: " + msg + "\n")
    raise SystemExit(1)


def path(n):
    return sparse.diags([np.ones(n - 1), np.ones(n - 1)], [-1, 1], format="csr")


def poly_of(T, coeffs):
    """coeffs[k] multiplies T^k (k=0..D); Horner in sparse arithmetic."""
    n = T.shape[0]
    M = sparse.csr_matrix((n, n))
    for c in reversed(coeffs):
        M = M @ T + c * sparse.identity(n, format="csr")
    return M


def block_hamiltonian(L, coeffs):
    """Principal compression of f(T_2L) to the first L sites.  Paths of length <= D starting and
    ending in the block never pass site L+D, so the compression of f(T_{L+D}) is exact."""
    D = len(coeffs) - 1
    n = min(2 * L, L + D)
    M = poly_of(path(n), coeffs).toarray()
    return M[:L, :L]


def block_correlation(L):
    """Positive-energy Fermi projector of T_2L restricted to the first L sites (App. A formula)."""
    N = 2 * L
    i = np.arange(1, L + 1)

    def Dn(r):
        r = np.asarray(r, dtype=float)
        out = np.empty_like(r)
        zero = r == 0
        out[zero] = 0.5
        rr = r[~zero]
        out[~zero] = np.sin(np.pi * rr / 2) / (2 * (N + 1) * np.sin(np.pi * rr / (2 * (N + 1))))
        return out

    I, J = np.meshgrid(i, i, indexing="ij")
    return Dn(I - J) - Dn(I + J)


def smallest_subset_sums(vals, M):
    """The M smallest subset sums of a nonnegative list, in nondecreasing order (with multiplicity)."""
    v = np.sort(np.asarray(vals, dtype=float))
    out = [0.0]
    if len(v) == 0:
        return np.array(out)
    heap = [(v[0], 0)]
    while heap and len(out) < M:
        s, k = heapq.heappop(heap)
        out.append(s)
        if k + 1 < len(v):
            heapq.heappush(heap, (s + v[k + 1], k + 1))
            heapq.heappush(heap, (s - v[k] + v[k + 1], k + 1))
    return np.array(out)


def unrestricted_Q(x, a, M):
    """x: entanglement energies |log((1-nu)/nu)|, a: one-body costs |eps|.  Returns Q and the
    certified bound on the truncation error."""
    xs = np.minimum(x, 700.0)
    logZ = float(np.sum(np.log1p(np.exp(-xs))))
    sx = smallest_subset_sums(xs, M)
    se = smallest_subset_sums(a, M)
    m = min(len(sx), len(se))
    p = np.exp(-sx[:m] - logZ)
    Q = float(np.dot(p, se[:m]))
    missing = max(0.0, 1.0 - float(np.sum(p)))
    # omitted ranks carry energy at most the full bandwidth sum(a)
    return Q, missing, missing * float(np.sum(a))


def analyse(L, coeffs, M):
    B = block_hamiltonian(L, coeffs)
    eps = eigvalsh(B)
    C = block_correlation(L)
    nu = np.clip(eigvalsh(C), 1e-300, 1 - 1e-16)
    with np.errstate(divide="ignore"):
        x = np.abs(np.log(nu) - np.log1p(-nu))
    h = -(nu * np.log(nu) + (1 - nu) * np.log1p(-nu))
    S = float(np.sum(np.nan_to_num(h)))
    a = np.abs(eps)
    Q, missing, err = unrestricted_Q(x, a, M)
    # Gaussian number-conserving: largest occupations on most negative energies
    occ = np.sort(1.0 - nu)[::-1]  # filling negative modes: occupation 1-nu_pos
    QG = float(np.dot(np.sort(eps), occ) - np.sum(eps[eps < 0]))
    # charge distribution of the block (Poisson binomial) and forced-cost floor
    neg = np.sort(-eps[eps < -1e-12])
    pos = np.sort(eps[eps > 1e-12])
    nstar = len(neg)
    nzero = L - len(neg) - len(pos)
    dist = np.zeros(1)
    dist[0] = 1.0
    for q in occ:
        dist = np.convolve(dist, [1 - q, q])
    cneg = np.concatenate([[0.0], np.cumsum(neg)])
    cpos = np.concatenate([[0.0], np.cumsum(pos)])
    floor = 0.0
    for n, pn in enumerate(dist):
        if pn < 1e-300:
            continue
        k = nstar - n
        if k > 0:
            floor += pn * cneg[min(k, len(neg))]
        j = n - nstar - nzero
        if j > 0:
            floor += pn * cpos[min(j, len(pos))]
    if not (floor <= QG + 1e-12):
        die("floor exceeds QG at L=%d" % L)
    if not (Q <= QG + err + 1e-12):
        die("Q exceeds QG at L=%d" % L)
    return dict(S=S, Q=Q, QG=QG, floor=floor, missing=missing, err=err)


def constants():
    eta = lambda s: (1 - mp.power(2, 1 - s)) * mp.zeta(s)
    out = {}
    for z in (1, 3):
        Az = mp.gamma(1 + mp.mpf(1) / z) * eta(1 + mp.mpf(1) / z) / mp.pi
        Dz = mp.power(z, z) / mp.power(z + 1, z + 1) / mp.power(Az, z)
        Gz = mp.gamma(z + 1) * eta(z + 1) / mp.power(mp.pi, z + 2)
        out[z] = (Dz, Gz)
    a = mp.mpf(3) / 2
    kap = a * (4 - a * a)
    A = 2 * mp.gamma(mp.mpf(3) / 2) * eta(mp.mpf(3) / 2) / (mp.pi * mp.sqrt(kap))
    out["three"] = (4 / (27 * A * A), a * a / (4 * mp.pi), a * a / (24 * mp.pi))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="16,64,256,1024,4096")
    ap.add_argument("--ranks", type=int, default=65536)
    args = ap.parse_args()
    sizes = [int(s) for s in args.sizes.split(",")]
    K = constants()
    D1, G1 = K[1]
    D3, G3 = K[3]
    c1, c2, c3 = K["three"]
    print("model,L,S,Q,QG,charge_floor,Q_scaled,QG_scaled,floor_scaled,QG_over_Q,"
          "limit_Q_scaled,limit_QG_scaled,limit_floor_scaled,limit_QG_over_Q,omitted_mass,omitted_energy_bound")
    for L in sizes:
        d = math.log(L)
        # z = 1 and z = 3 odd-power chains (same state)
        for z, Dz, Gz in ((1, D1, G1), (3, D3, G3)):
            coeffs = [0.0] * z + [1.0]
            r = analyse(L, coeffs, args.ranks)
            S = r["S"]
            qs = r["Q"] * L ** z / S ** (z + 1)
            gs = r["QG"] * L ** z / d ** (z + 1)
            print("z=%d,%d,%.10f,%.12e,%.12e,,%.8f,%.8f,,%.8f,%.8f,%.8f,,%.8f,%.2e,%.2e" % (
                z, L, S, r["Q"], r["QG"], qs, gs, r["QG"] / r["Q"], float(Dz), float(Gz),
                float(Gz / (Dz / 6 ** (z + 1))), r["missing"], r["err"]))
            sys.stdout.flush()
        # range-three chain f(x) = x (x - 3/2)^2 = x^3 - 3x^2 + 9x/4
        r = analyse(L, [0.0, 2.25, -3.0, 1.0], args.ranks)
        S = r["S"]
        print("three,%d,%.10f,%.12e,%.12e,%.12e,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,,%.2e,%.2e" % (
            L, S, r["Q"], r["QG"], r["floor"], r["Q"] * L ** 2 / S ** 3, r["QG"] * L / d ** 2,
            r["floor"] * L / d, r["QG"] / r["Q"], float(c1), float(c3), float(c2), r["missing"], r["err"]))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
