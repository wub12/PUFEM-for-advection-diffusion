import matplotlib.pyplot as plt
import numpy as np


#B = np.loadtxt('result.txt', delimiter=',', usecols=range(320));
B = np.loadtxt('error.txt', delimiter=',', usecols=range(320));

plt.figure(2)
plt.imshow(B, interpolation='none',  extent=[0, 1, 0, 1])
plt.grid(True)
#plt.imshow(B, interpolation='none',  extent=[0, 1, 0, 1])
#plt.grid(True)

plt.show()