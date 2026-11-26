"""Plot the two zero-energy points of the cubic hopping polynomial."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
a=1.5
k=np.linspace(0.0, 0.60, 1000)
x=2*np.cos(np.pi*k)
energy=x*(x-a)**2
touching=np.arccos(a/2)/np.pi
fig, ax=plt.subplots(figsize=(6,2.55), layout="constrained")
ax.axhline(0, color="#8d827e", lw=.8)
ax.axvspan(.5,.6,color="#eee4d3")
ax.plot(k,energy,lw=2.5,color="#70405c")
ax.plot([touching,.5],[0,0],"o",color="#70405c",ms=4)
ax.annotate("Quadratic particle modes",xy=(touching,.04),xytext=(.115,.70),
            arrowprops={"arrowstyle":"->","color":"#643953"},fontsize=9)
ax.annotate("Linear crossing",xy=(.5,.015),xytext=(.413,.70),
            arrowprops={"arrowstyle":"->","color":"#643953"},fontsize=9)
ax.text(.553,-.48,"Filled",color="#6f5635",fontsize=9,ha="center")
ax.spines[["top","right"]].set_visible(False)
ax.set(xlim=(0,.6),ylim=(-.64,.9),xlabel=r"Wave number $k/\pi$",ylabel="One-particle energy")
fig.savefig(Path(__file__).resolve().parent.parent / "figures" / "dispersion.pdf")
