"""Numerical check of the pressure constant in the general dispersion theorem (Theorem 5).

For a sign-preserving polynomial f, the one-body costs of the open path are |f(2 cos k_j)|,
k_j = pi j/(L+1).  Theorem 5 states that, at inverse temperature beta = b (L/u)^m with
1 << u << L, log Z / u -> A b^(-1/m), where m is the largest order of a zero of f in (-2,2) and
    A = (2/pi) Gamma(1+1/m) eta(1+1/m) * sum over zeros of order m of kappa^(-1/m),
    kappa = |c| (4 - a^2)^(m/2)   for f(x) ~ c (x-a)^m near the zero a.
The theorem itself uses u = S_L ~ (log L)/6; this script uses u = sqrt(L) only so that the
Riemann sum is resolved at a modest L.  It also checks that the Legendre transform gives
D(m,A) = m^m/((m+1)^(m+1) A^m) and reproduces the constants quoted in the manuscript
D_1, D_3, the range-three coefficient, R_3, R_5, the large-z ratio and the pair ratio.
Usage: python3 code/general_dispersion_check.py > results/general_dispersion.csv
"""
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 30


def require(condition, message):
    """Explicit check that survives python -O (no assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED: " + message)


def eta(s):
    return (1 - mp.power(2, 1 - s)) * mp.zeta(s)


def A_of(kappas, m):
    return 2 * mp.gamma(1 + mp.mpf(1) / m) * eta(1 + mp.mpf(1) / m) / mp.pi * sum(
        mp.power(k, -mp.mpf(1) / m) for k in kappas)


def D_of(m, A):
    return mp.power(m, m) / mp.power(m + 1, m + 1) / mp.power(A, m)


CASES = [
    # name, f, m, list of kappa for the zeros of maximal order, tolerance at L = 1e7
    ("x^3", lambda x: x ** 3, 3, [mp.mpf(8)], 1e-4),
    ("x(x-3/2)^2", lambda x: x * (x - 1.5) ** 2, 2, [mp.mpf("1.5") * (4 - mp.mpf("2.25"))], 1e-4),
    ("x(x-3/2)^2(x+1/2)^2", lambda x: x * (x - 1.5) ** 2 * (x + 0.5) ** 2, 2,
     [6 * (4 - mp.mpf("2.25")), 2 * (4 - mp.mpf("0.25"))], 3e-3),
    ("x^3(x-3/2)^2", lambda x: x ** 3 * (x - 1.5) ** 2, 3, [mp.mpf("2.25") * 8], 3e-2),
]


def main():
    print("case,m,L,b,logZ_over_u,predicted,ratio")
    for name, f, m, kappas, tol in CASES:
        A = A_of(kappas, m)
        for L in (10 ** 5, 10 ** 6, 10 ** 7):
            k = np.pi * np.arange(1, L + 1) / (L + 1)
            e = np.abs(f(2 * np.cos(k)))
            u = np.sqrt(L)
            for b in (0.5, 1.0, 2.0):
                beta = b * (L / u) ** m
                lz = float(np.sum(np.log1p(np.exp(-beta * e))))
                pred = float(A * mp.power(b, -mp.mpf(1) / m))
                ratio = lz / u / pred
                print("%s,%d,%d,%g,%.8f,%.8f,%.6f" % (name, m, L, b, lz / u, pred, ratio))
                if L == 10 ** 7:
                    require(abs(ratio - 1) < tol, "%s pressure ratio %.6f" % (name, ratio))
    # Legendre transform: inf_b (b e + A b^(-1/m)) = 1 exactly at e = D(m, A)
    for m in (1, 2, 3, 5):
        A = mp.mpf("0.37")
        e = D_of(m, A)
        val = mp.power(A / (m * e), mp.mpf(m) / (m + 1))  # stationary point of b e + A b^(-1/m)
        sigma = val * e + A * mp.power(val, -mp.mpf(1) / m)
        require(abs(sigma - 1) < mp.mpf(10) ** -25, "Legendre identity m=%d" % m)
    # constants quoted in the manuscript
    D1 = D_of(1, A_of([mp.mpf(2)], 1))
    D3 = D_of(3, A_of([mp.mpf(8)], 3))
    a = mp.mpf(3) / 2
    three = D_of(2, A_of([a * (4 - a * a)], 2))
    require(abs(D1 - 3 / mp.pi) < mp.mpf(10) ** -25, "D_1 = 3/pi")
    require(abs(D3 - mp.mpf("11.2023136908")) < mp.mpf(10) ** -9, "D_3")
    require(abs(three - mp.mpf("2.086823")) < mp.mpf(10) ** -6, "range-three coefficient")
    # one-body (Gaussian) to unrestricted limit ratios R_z = G_z 6^(z+1) / D_z (Theorem 11)
    def R(z):
        Dz = D_of(z, A_of([mp.power(2, z)], z))
        Gz = mp.gamma(z + 1) * eta(z + 1) / mp.power(mp.pi, z + 2)
        return Gz * mp.power(6, z + 1) / Dz
    R1, R3, R5 = R(1), R(3), R(5)
    require(abs(R1 - 1) < mp.mpf(10) ** -25, "R_1 = 1")
    require(abs(R3 - mp.mpf("2.14815")) < mp.mpf("5e-6"), "R_3")
    require(abs(R5 - mp.mpf("11.55617")) < mp.mpf("5e-6"), "R_5")
    zbig = 2001
    large = mp.power(R(zbig), mp.mpf(1) / zbig) / zbig
    limit = 6 * mp.log(2) / (mp.e * mp.pi ** 2)
    require(abs(large / limit - 1) < mp.mpf("0.01"), "large-z ratio approaches 6 log2/(e pi^2)")
    # pair-creating Gaussian control in the range-three chain (Section 7): ratio to Q_L
    kap = a * (4 - a * a)
    Apair = A_of([kap], 2)
    pair = kap * mp.gamma(3) * eta(3) / (4 * mp.pi ** 4)
    unres = D_of(2, Apair) / mp.mpf(6) ** 3  # S_L ~ d/6
    ratio_pair = pair / unres
    require(abs(ratio_pair - mp.mpf("1.25734")) < mp.mpf("5e-6"), "pair ratio")
    sys.stderr.write("general_dispersion_check: all checks passed; D_1=%s D_3=%s three=%s "
                     "R_3=%s R_5=%s pair=%s large-z(%d)=%s limit=%s\n"
                     % (mp.nstr(D1, 12), mp.nstr(D3, 12), mp.nstr(three, 10), mp.nstr(R3, 8),
                        mp.nstr(R5, 9), mp.nstr(ratio_pair, 8), zbig, mp.nstr(large, 6),
                        mp.nstr(limit, 6)))


if __name__ == "__main__":
    main()
