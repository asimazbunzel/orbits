
import numpy as np
import matplotlib.pyplot as plt

fname = "orbit_history.data"
orbit = np.loadtxt(fname, skiprows=True, unpack=True)

ax = plt.gca()
ax.set_aspect("equal", "datalim")

ax.set_xlabel('X [R$_\\odot$]')
ax.set_ylabel('Y [R$_\\odot$]')

ax.plot(orbit[3], orbit[4])

ax.scatter([0], [0], s=250, marker="*", edgecolor="black", facecolor="yellow")

plt.show()
