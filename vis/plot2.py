
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

data  = np.loadtxt("../data/f.dat", delimiter=",")
data2 = np.loadtxt("../data/g.dat", delimiter=",")

norm = Normalize(vmin=-0.5, vmax=0.)

fig, axes = plt.subplots(
    2, 1,
    sharex=True,
    sharey=True,
    figsize=(6, 6),
    constrained_layout=True   # ← key change
)

im1 = axes[0].imshow(data, origin="lower", aspect="auto", norm=norm)
axes[0].set_title("f.dat")

im2 = axes[1].imshow(data2, origin="lower", aspect="auto", norm=norm)
axes[1].set_title("g.dat")
axes[1].set_xlabel("Column index")

fig.colorbar(im1, ax=axes, label="Value (0–1)")

plt.show()

