#!/usr/bin/env python
#
import os


import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from mpl_toolkits.axes_grid.inset_locator import inset_axes

dirbase1 = "0823_0"
dirbase2 = "0823_1"

input1 = os.listdir(dirbase1 + "/" + dirbase1 + "_ABUND")

input2 = os.listdir(dirbase2 + "/" + dirbase2 + "_ABUND")

try:
    input1.remove('.DS_Store')
    input2.remove('.DS_Store')
except Exception, e:
    pass


properties = np.array(np.loadtxt("properties.txt", dtype={'names': ('label', 'time', 'T(GK)', 'rho(cgs)'),
                                                          'formats': (np.int, np.float, np.float, np.float)}).tolist()).transpose()
# print properties
time = []
temperature = []
density = []

# for fname in finput:
for ii in xrange(0, len(input1)):
    fname1 = input1[ii]
    fname2 = input2[ii]
    time.append(properties[1][ii])
    temperature.append(properties[2][ii])
    density.append(properties[3][ii])
    # print fname
    c = np.loadtxt("./dot_abund/" + fname, dtype={'names': ('a', 'abund', 'mass frac', 'norm abund'),
                                                  'formats': (np.int, np.float, np.float, np.float)})
    c = np.array(c.tolist()).transpose()
    plt.clf()
    plt.yscale('log', nonposy='clip')
    plt.ylabel("Abundance per nucleon")
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
    plt.savefig("dot_abund_png/" + fname + ".png")
