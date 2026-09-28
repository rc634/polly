import numpy as np
import matplotlib.pyplot as plt

data_dir = "../data"

xL, xU = 0., 40.
yL, yU = 0., 40.
extent = [-xU, xU, -yU, yU]
physical_aspect = 1.0

def reflect_full(arr):
    arr_x = np.hstack([np.fliplr(arr), arr])
    return np.vstack([np.flipud(arr_x), arr_x])

files = {
    "source-W":   f"{data_dir}/source-W.dat",
    "source-rho": f"{data_dir}/source-psi.dat",
    "W":          f"{data_dir}/W.dat",
    "psi":        f"{data_dir}/psi.dat",
}

arrays = {}
for name, path in files.items():
    arr = np.loadtxt(path, delimiter=",")
    arrays[name] = arr
    print(f"  {name:12s}  shape = {arr.shape[0]} x {arr.shape[1]}  (ny x nx)")

# row 0: W source | W solution
# row 1: rho source | psi solution
datasets = [
    (arrays["source-W"],   "W source (rotation)"),
    (arrays["W"],          "W solution"),
    (arrays["source-rho"], "density source (rho)"),
    (arrays["psi"],        "psi (conformal factor)"),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)

plot_order = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]

for ax, (data, title) in zip(plot_order, datasets):
    im = ax.imshow(
        reflect_full(data),
        origin="lower",
        aspect=physical_aspect,
        extent=extent,
        cmap="rainbow",
    )
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim([-4, 4])
    ax.set_ylim([-4, 4])
    fig.colorbar(im, ax=ax)

fig.text(0.25, 1.01, "sources",   ha="center", fontsize=11, color="grey", transform=axes[0,0].transAxes)
fig.text(0.25, 1.01, "solutions", ha="center", fontsize=11, color="grey", transform=axes[0,1].transAxes)

plt.show()
