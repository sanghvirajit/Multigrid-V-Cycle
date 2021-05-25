import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

Y1 = []

X = np.linspace(0, 1, 129)

#PLOT Y
y2 = pd.read_csv("w.txt", header=None)

for i in range(0, len(X)):
	x = X[i]
	y = x * (1-x)**2
	Y1.append(y)

fig = plt.gcf()

plt.plot(X, Y1, "k+")
plt.plot(X, Y1, label="y analytical")

plt.plot(X, y2, "kx")
plt.plot(X, y2, label="y exp")

plt.xlabel("Space Grid")
plt.ylabel("y")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical.png")
plt.show()


