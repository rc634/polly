import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Load data
# data = np.loadtxt("../data/psi.dat", delimiter=",")
# data2 = np.loadtxt("../data/dataPsiNeumann.csv", delimiter=",")
data = np.loadtxt("../data/W.dat", delimiter=",")
data2 = np.loadtxt("../data/dataWNeumann2.csv", delimiter=",")
# data2 = np.transpose(data2)

# Optional normalization
norm = Normalize(vmin=-1, vmax=0.)

fig, axes = plt.subplots(
    2, 1,
    sharex=True,
    sharey=True,
    figsize=(6, 6),
    constrained_layout=True
)

# First image (e.g. plasma)
im1 = axes[0].imshow(
    np.log10(np.abs(data)),
    origin="lower",
    aspect="auto",
    cmap="rainbow"      # <-- colormap 1
    # norm=norm
)
axes[0].set_title("log10|psi-1|.dat")
cbar1 = fig.colorbar(im1, ax=axes[0], label="Value")

# Second image (e.g. viridis)
im2 = axes[1].imshow(
    np.log10(np.abs(data2)),
    origin="lower",
    aspect="auto",
    cmap="rainbow",     # <-- colormap 2
    # norm=norm
)
axes[1].set_title("log10|W|.dat")
axes[1].set_xlabel("Column index")
cbar2 = fig.colorbar(im2, ax=axes[1], label="Value")

plt.show()
