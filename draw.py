import numpy as np
import matplotlib.pyplot as plt
import csv

x = []
y = []

x, y = np.loadtxt('out.txt', delimiter=',', unpack=True)

maxX = max(x)
minX = min(x)
maxY = max(y)
minY = min(y)

plt.axis([-1000, 1000, -1000, 1000])
plt.grid(True)
plt.scatter(x, y, label='collision', s=50)
plt.title("Rutherford Results")

plt.legend()
plt.show()