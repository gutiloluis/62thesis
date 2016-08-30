import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt("one_cell.dat")

plt.plot(A[:,0],A[:,1])
plt.plot(A[:,0],A[:,3])
plt.show()
