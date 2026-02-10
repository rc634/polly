import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

data  = np.loadtxt("../data/W.dat", delimiter=",")
data2 = np.loadtxt("../data/psi.dat", delimiter=",")

fig, axes = plt.subplots(
    2, 1,
    sharex=True,
    sharey=True,
    figsize=(6, 6),
    constrained_layout=True
)

# ---- W plot ----
im1 = axes[0].imshow(data, origin="lower", aspect="auto")
axes[0].set_title("W.dat")
cbar1 = fig.colorbar(im1, ax=axes[0])
cbar1.set_label("W value")

# ---- psi plot ----
im2 = axes[1].imshow(data2, origin="lower", aspect="auto")
axes[1].set_title("psi.dat")
axes[1].set_xlabel("Column index")
cbar2 = fig.colorbar(im2, ax=axes[1])
cbar2.set_label("psi value")

# Add isocontour psi = 2
axes[1].contour(
    data2,
    levels=[2.0],
    colors="black",
    linewidths=1.5,
    origin="lower"
)

plt.show()
