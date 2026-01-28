import numpy as np
import matplotlib.pyplot as plt

# load the data (comma-separated)
data = np.loadtxt("../data/f.dat", delimiter=",")
data2 = np.loadtxt("../data/g.dat", delimiter=",")

plt.imshow(data-data2, origin="lower", aspect="auto")
plt.colorbar(label="Value")
plt.xlabel("Column index")
plt.ylabel("Row index")

plt.show()
