import numpy as np
import matplotlib.pyplot as plt

# Physical domain boundaries (match params.hpp)
xL, xU = 0., 10.
yL, yU = 0., 10.
extent = [xL, xU, yL, yU]
physical_aspect = (yU - yL) / (xU - xL)

# Switch: True = log plots for psi and W, False = raw values
logarithmic_metric = False

# Load data
files = {
    "rho-raw":    "../data/rho-raw.dat",
    "rho-eff":    "../data/rho-eff.dat",
    "source-psi": "../data/source-psi.dat",
    "psi":        "../data/psi.dat",
    "v-raw":      "../data/v-raw.dat",
    "v-eff":      "../data/v-eff.dat",
    "source-W":   "../data/source-W.dat",
    "W":          "../data/W.dat",
}

arrays = {}
for name, path in files.items():
    arr = np.loadtxt(path, delimiter=",")
    arrays[name] = arr
    print(f"  {name:12s}  shape = {arr.shape[0]} x {arr.shape[1]}  (ny x nx)")

if logarithmic_metric:
    psi_data  = np.log10(np.abs(arrays["psi"] - 1))
    psi_title = "log10|psi - 1|"
    W_data    = np.log10(np.abs(arrays["W"] + 1e-25))
    W_title   = "log10|W|"
else:
    psi_data  = arrays["psi"]
    psi_title = "psi"
    W_data    = arrays["W"]
    W_title   = "W"

datasets = [
    # row 1: psi / rho
    (arrays["rho-raw"],    "rho raw"),
    (arrays["rho-eff"],    "rho eff"),
    (arrays["source-psi"], "source psi"),
    (psi_data,             psi_title),
    # row 2: W / v
    (arrays["v-raw"],      "v raw"),
    (arrays["v-eff"],      "v eff"),
    (arrays["source-W"],   "source W"),
    (W_data,               W_title),
]

fig, axes = plt.subplots(
    2, 4,
    figsize=(18, 8),
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