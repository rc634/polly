import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# ---- Set physical domain size here ----
H = 10.0   # height (y max)
W = 10.0    # width  (x max)

data  = np.loadtxt("../data/W.dat", delimiter=",")
data2 = np.loadtxt("../data/psi.dat", delimiter=",")

fig, axes = plt.subplots(
    2, 1,
    sharex=True,
    sharey=True,
    figsize=(6, 6),
    constrained_layout=True
)

# Define extent: [xmin, xmax, ymin, ymax]
extent = [0, W, 0, H]

# ---- W plot ----
im1 = axes[0].imshow(
    np.log(data*data),
    origin="lower",
    aspect="auto",
    extent=extent
)
axes[0].set_title("W.dat")
axes[0].set_ylabel("z")
# axes[0].set_xlabel("cylindrical radius")
cbar1 = fig.colorbar(im1, ax=axes[0])
cbar1.set_label(r"$\ln(W^2)$")

# ---- psi plot ----
im2 = axes[1].imshow(
    np.log(data2-1.),
    origin="lower",
    aspect="auto",
    extent=extent
)
axes[1].set_title("psi.dat")
axes[1].set_ylabel("z")
axes[1].set_xlabel("cylindrical radius")

cbar2 = fig.colorbar(im2, ax=axes[1])
cbar2.set_label(r"$\ln(\psi - 1)$")

# # Add isocontour psi = 2
# axes[1].contour(
#     np.linspace(0, W, data2.shape[1]),
#     np.linspace(0, H, data2.shape[0]),
#     data2,
#     levels=[2.0],
#     colors="black",
#     linewidths=1.5
# )

for ax in axes:
    ax.set_aspect('equal')

plt.show()