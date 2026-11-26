"""Evaluate the cubic crossover for the unit-hopping open chain."""
from pathlib import Path
import csv
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit, gamma, zeta
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def require(condition, message):
    """Explicit check that survives python -O (no assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED: " + message)


OUT = Path(__file__).resolve().parent.parent  # package root

def crossover(z, lam):
    # Scale x so that b times the dispersion is order one in the integral.
    def thermodynamics(log_b):
        b = np.exp(log_b)
        x0 = min(b**(-1/z), 1/(b*lam)) if lam else b**(-1/z)
        def integrand(y, moment):
            w = b*((x0*y)**z + lam*x0*y)
            return w*expit(-w) if moment else np.logaddexp(0, -w)
        p, mean = (x0/np.pi*quad(integrand, 0, np.inf, args=(m,), epsabs=3e-12)[0]
                   for m in (False, True))
        return p + mean, mean/b
    root = brentq(lambda h: thermodynamics(h)[0] - 1, -25, 20, xtol=2e-13)
    return np.exp(root), thermodynamics(root)[1]

def main():
    z = 3
    a = gamma(1+1/z)*(1-2**(-1/z))*zeta(1+1/z)/np.pi
    d = z**z/((z+1)**(z+1)*a**z)
    require(abs(crossover(z, 0)[1]/d-1)<1e-11, "crossover.py:32")
    ls = np.logspace(-2, 4, 73)
    vals = np.array([crossover(z, x) for x in ls])
    require(np.all(np.diff(vals[:,1])>0), "crossover.py:35")
    require(abs(vals[-1,1]/(ls[-1]*3/np.pi)-1)<0.003, "crossover.py:36")
    with (OUT/"results"/"crossover.csv").open("w", newline="") as handle:
        writer=csv.writer(handle)
        writer.writerow(["lambda","b","scaled_passive_energy"])
        writer.writerow([0,*crossover(z,0)])
        writer.writerows(zip(ls, *vals.T))
    fig, ax=plt.subplots(figsize=(6.8,3.1), layout="constrained")
    ax.loglog(ls, vals[:,1], lw=2.6, color="#71395b", label="Cubic chain with linear hopping")
    ax.axhline(d, ls="--", lw=1.3, color="#aa642d", label="Pure cubic limit")
    ax.loglog(ls, 3*ls/np.pi, ls=":", lw=1.8, color="#247973", label="Linear limit")
    ax.set(xlabel=r"Scaled linear hopping  $\lambda=a(L/S)^2$", ylabel=r"Scaled residual energy  $Q L^3/S^4$",
           ylim=(1,2e4))
    ax.legend(frameon=False, loc="upper left", fontsize=8.5)
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", which="major", color="#e3ddda", linewidth=.6)
    fig.savefig(OUT/"figures"/"crossover.pdf")
    print(f"Pure cubic coefficient {d:.12f}")
    print("Cubic crossover at lambda=1:",crossover(3,1))
    print("Linear limit ratio at lambda=10000:",vals[-1,1]/(ls[-1]*3/np.pi))
if __name__=="__main__":
    main()
