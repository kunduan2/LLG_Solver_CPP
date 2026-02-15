import numpy as np
import matplotlib.pyplot as plt

# Load data: columns => t, mx, my, mz
data = np.loadtxt("test.dat")

t = data[:,0]
mx = data[:,1]
my = data[:,2]
mz = data[:,3]

plt.figure(figsize=(8,5))

plt.plot(t, mx, linewidth=2, label="mx")
plt.plot(t, my, linewidth=2, label="my")
plt.plot(t, mz, linewidth=2, label="mz")

plt.xlabel("Time")
plt.ylabel("Magnetization")
plt.title("LLG Magnetization Dynamics")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
