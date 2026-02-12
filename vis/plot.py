import numpy as np
import matplotlib.pyplot as plt

# load the data (comma-separated)
data = np.loadtxt("../data/psi.dat", delimiter=",")

plt.imshow(np.log(np.abs(data-1)), origin="lower", aspect="auto")
plt.colorbar(label="Value")
plt.xlabel("Column index")
plt.ylabel("Row index")

plt.show()
