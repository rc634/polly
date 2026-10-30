import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 15,
    "legend.fontsize": 12,
    "lines.linewidth": 2.0,
    "mathtext.fontset": "cm",
})

def f(x, s):
    return np.exp(s**2 * x**2 / (x**2 - 1))

x = np.linspace(-1, 1, 2000)
s_values = [0., 1./3., 1., np.sqrt(3.), 3.]

fig, ax = plt.subplots(figsize=(10, 4.5))


ax.plot(x, f(x, s_values[0]), label=r"$\sigma = 0$")
ax.plot(x, f(x, s_values[1]), label=r"$\sigma = 1/3$")
ax.plot(x, f(x, s_values[2]), label=r"$\sigma = 1$")
ax.plot(x, f(x, s_values[3]), label=r"$\sigma = \sqrt{3}$")
ax.plot(x, f(x, s_values[4]), label=r"$\sigma = 3$")

ax.set_xlim(-1.1, 1.1)
ax.set_ylim(0, 1.1)
ax.set_xlabel(r"$\xi$")
ax.set_ylabel(r"$P(\xi;\sigma)$")
ax.set_title(r"Bump Function")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)

fig.savefig("bump_plot.pdf")
fig.savefig("bump_plot.png")
plt.show()