from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

a0001 = np.array(np.loadtxt('mg25_0.01.dat', dtype={'names': ('label', 'abund'),
                                                    'formats': (np.int, np.float)}).tolist()).transpose()
a01 = np.array(np.loadtxt('mg25_0.1.dat', dtype={'names': ('label', 'abund'),
                                                 'formats': (np.int, np.float)}).tolist()).transpose()
a1 = np.array(np.loadtxt('mg25_1.dat', dtype={'names': ('label', 'abund'),
                                              'formats': (np.int, np.float)}).tolist()).transpose()
a10 = np.array(np.loadtxt('mg25_10.dat', dtype={'names': ('label', 'abund'),
                                                'formats': (np.int, np.float)}).tolist()).transpose()
a100 = np.array(np.loadtxt('mg25_100.dat', dtype={'names': ('label', 'abund'),
                                                  'formats': (np.int, np.float)}).tolist()).transpose()
a1000 = np.array(np.loadtxt('mg25_1000.dat', dtype={'names': ('label', 'abund'),
                                                    'formats': (np.int, np.float)}).tolist()).transpose()

x = np.array([0.01, 0.1, 1, 10, 100, 1000])
y = np.array([a0001[1][len(a0001[1]) - 1], a01[1][len(a01[1]) - 1], a1[1][len(a1[1]) - 1], a10[1][len(a10[1]) - 1],
              a100[1][len(a100[1]) - 1], a1000[1][len(a1000[1]) - 1]])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.xlabel(r'Factor')
plt.ylabel(r"Mass fraction of {}^{25}Mg")
# print x, y
plt.plot(x, y, 'o-')
# plt.yscale('log')
plt.xscale('log')

plt.show()
