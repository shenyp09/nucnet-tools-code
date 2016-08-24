#!/usr/bin/env python
#
import os
import re

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from mpl_toolkits.axes_grid.inset_locator import inset_axes

if len(argv) != 3:
    print "\nUsage: %s propertyfile abundfile\n" % argv[0]
    exit(1)

script, propertyfile, abundfile = argv
label = int(re.findall(r"\d+", abundfile)[0])

properties = np.array(np.loadtxt(propertyfile, dtype={'names': ('label', 'time', 'T(GK)', 'rho(cgs)'),
                                                      'formats': (np.int, np.float, np.float, np.float)}).tolist()).transpose()
time = properties[1][:label]
temperature = properties[2][:label]
density = properties[3][:label]

# print fname
c = np.loadtxt(abundfile, dtype={'names': ('a', 'abund'),
                                 'formats': (np.int, np.float)})
c = np.array(c.tolist()).transpose()
plt.clf()
plt.yscale('log', nonposy='clip')
plt.ylabel("Abundance")
plt.xlabel("A")
axes = plt.gca()
axes.set_ylim([1e-10, 1])
plt.plot(c[0], c[1], ls='-', color='black',
         marker='.', mec='blue', mew=2., ms=10.)
# this is an inset axes over the main axes
inset = inset_axes(axes,
                   width="40%",  # width = 30% of parent_bbox
                   height=1,  # height : 1 inch
                   loc=1)
# n, bins, patches = plt.hist(s, 400, normed=1)
plt.ylabel('T(GK)')
plt.xlabel('time(s)')
# print properties[1]
plt.plot(time, temperature)
plt.yscale('log', nonposy='clip')
plt.xscale('log', nonposy='clip')
ax = plt.gca()
ax.set_xlim([1e-5, 1e16])
ax.set_ylim([1e-7, 1])
plt.xticks(rotation=45)
ax2 = plt.twinx()  # this is the important function
ax2.plot(time, density)
ax2.set_yscale('log')
ax2.set_xlim([1e-5, 1e16])
ax2.set_ylim([1e-12, 1e4])
ax2.set_ylabel('$\\rho$(cgs)')
# n, bins, patches = plt.hist(s, 400, normed=1)
# plt.xticks(100)
# plt.yticks([])
plt.savefig("dot_abund_png/abund_%05d.dat.png" % label)
