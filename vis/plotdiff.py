import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# load the data (comma-separated)
data = np.loadtxt("../data/f.dat", delimiter=",")
data2 = np.loadtxt("../data/g.dat", delimiter=",")

# norm = Normalize(vmin=-5, vmax=0.)

# plt.imshow((((data-data2)/data2)), origin="lower", aspect="auto")#, norm=norm)
plt.imshow(np.log10(np.abs((data-data2)/data2)), origin="lower", aspect="auto")#, norm=norm)
# plt.imshow(data-data2, origin="lower", aspect="auto", norm=norm)
plt.colorbar(label="Log Fractional error")
plt.xlabel("gridpoints")
plt.ylabel("gridpoints")

plt.show()
