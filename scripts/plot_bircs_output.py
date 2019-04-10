#!/usr/bin/env python3
from mpl_toolkits.mplot3d import Axes3D 
import sys
import numpy as np
from matplotlib import pyplot as plt

#r, theta, psi
rtp = []

fp = open(sys.argv[1])
for ln in fp:
    bits = ln.split(',')
    rtp.append((float(bits[0].strip()), float(bits[2].strip()), float(bits[1].strip())))

fp.close()

    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs = []
ys = []
zs = []

for s in rtp:
    xs.append(s[0] * np.sin(s[1]) * np.cos(s[2]))
    ys.append(s[0] * np.sin(s[1]) * np.sin(s[2]))
    zs.append(s[0] * np.cos(s[1]))

ax.scatter(xs, ys, zs)

ax.set_xlabel('X RCS')
ax.set_ylabel('Y RCS')
ax.set_zlabel('Z RCS')

plt.show()
