import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# load the data (comma-separated)
data = np.loadtxt("../data/dataPsiDirichlet.csv", delimiter=",")
data2 = np.loadtxt("../data/dataPsiNeumann.csv", delimiter=",")
data = np.transpose(data)
data2 = np.transpose(data2)
# norm = Normalize(vmin=-5, vmax=0.)

# plt.imshow((((data-data2)/data2)), origin="lower", aspect="auto")#, norm=norm)
plt.imshow(np.log10(np.abs((data-data2)/data)), origin="lower", aspect="auto")#, norm=norm)
# plt.imshow(data-data2, origin="lower", aspect="auto", norm=norm)
plt.colorbar(label="Log Fractional error")
plt.xlabel("gridpoints")
plt.ylabel("gridpoints")

plt.show()
