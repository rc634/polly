import numpy as np
import matplotlib.pyplot as plt

# Physical domain boundaries (match params.hpp)
xL, xU = 0., 20.
yL, yU = 0., 20.
extent = [xL, xU, yL, yU]
physical_aspect = (yU - yL) / (xU - xL)

# Load data
files = {
    "source-psi":"../data/source-psi.dat",
    "psi":      "../data/psi.dat",
    "source-W": "../data/source-W.dat",
    "W":        "../data/W.dat",
}

arrays = {}
for name, path in files.items():
    arr = np.loadtxt(path, delimiter=",")
    arrays[name] = arr
    print(f"  {name:12s}  shape = {arr.shape[0]} x {arr.shape[1]}  (ny x nx)")

datasets = [
    (arrays["source-psi"],                        "source psi"),
    (np.log10(np.abs(arrays["psi"] - 1)),       "log10|psi - 1|"),
    (arrays["source-W"],                          "source W"),
    (np.log10(np.abs(arrays["W"]   + 1e-25)),   "log10|W|"),
]

fig, axes = plt.subplots(
    2, 2,
    figsize=(10, 8),
    constrained_layout=True
)

for ax, (data, title) in zip(axes.flat, datasets):
    im = ax.imshow(
        data,
        origin="lower",
        aspect=physical_aspect,
        extent=extent,
        cmap="rainbow",
    )
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax)

plt.show()