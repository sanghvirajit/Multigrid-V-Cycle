import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

Y1 = []

X = np.linspace(0, 1, 129)

for i in range(0, len(X)):
	x = X[i]
	y = x * (1-x)**2
	Y1.append(y)
#PLOT Y
y2 = pd.read_csv("w_tilda.txt", header=None)
y3 = pd.read_csv("w_tilda_vector.txt", header=None)
y4 = pd.read_csv("w_final.txt", header=None)

fig = plt.gcf()

#plt.plot(X, Y1, "k+")
plt.plot(X, Y1, label="w analytical")

#plt.plot(X, y2, "kx")
plt.plot(X, y2, label="w tilda")
plt.plot(X, y3, label="w tilda vector")
plt.plot(X, y4, label="w final")

plt.xlabel("Space Grid")
plt.ylabel("y")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical.png")
plt.show()


