import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import json

Y = []
F = []

X = np.linspace(0, 1, 65)

print(X)

for i in range(0, len(X)):
	x = X[i]
	y = x * (1-x)**2
	Y.append(y)
	
with open('Y.txt', 'w') as f:
	for item in Y:
        	f.write("%s\n" % item)
        
for i in range(0, len(X)):
	x = X[i]
	z = 4 - 6*x
	F.append(z)
	
with open('F.txt', 'w') as f:
	for item in F:
        	f.write("%s\n" % item)
