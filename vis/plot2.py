import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Load data
data = np.loadtxt("../data/psi.dat", delimiter=",")
data2 = np.loadtxt("../data/W.dat", delimiter=",")

# Normalization (optional)
norm = Normalize(vmin=-0.5, vmax=0.)

# Create subplots
fig, axes = plt.subplots(
    2, 1,
    sharex=True,
    sharey=True,
    figsize=(6, 6),
    constrained_layout=True
)

# First image
im1 = axes[0].imshow(np.log10(np.abs(data-1)), origin="lower", aspect="auto")#, norm=norm)
axes[0].set_title("log10|psi-1|.dat")
# Add colorbar for the first subplot
cbar1 = fig.colorbar(im1, ax=axes[0], label="Value (0–1)")

# Second image
im2 = axes[1].imshow(np.log10(np.abs(data2)), origin="lower", aspect="auto")#, norm=norm)
axes[1].set_title("log10|W|.dat")
axes[1].set_xlabel("Column index")
# Add colorbar for the second subplot
cbar2 = fig.colorbar(im2, ax=axes[1], label="Value (0–1)")

plt.show()

