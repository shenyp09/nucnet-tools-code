from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

a0001 = np.array(np.loadtxt('al26_0.01.dat', dtype={'names': ('label', 'abund', 'abund2'),
                                                    'formats': (np.int, np.float, np.float)}).tolist()).transpose()
a01 = np.array(np.loadtxt('al26_0.1.dat', dtype={'names': ('label', 'abund', 'abund2'),
                                                 'formats': (np.int, np.float, np.float)}).tolist()).transpose()
a1 = np.array(np.loadtxt('al26_1.dat', dtype={'names': ('label', 'abund', 'abund2'),
                                              'formats': (np.int, np.float, np.float)}).tolist()).transpose()
a10 = np.array(np.loadtxt('al26_10.dat', dtype={'names': ('label', 'abund', 'abund2'),
                                                'formats': (np.int, np.float, np.float)}).tolist()).transpose()
a100 = np.array(np.loadtxt('al26_100.dat', dtype={'names': ('label', 'abund', 'abund2'),
                                                  'formats': (np.int, np.float, np.float)}).tolist()).transpose()
a1000 = np.array(np.loadtxt('al26_1000.dat', dtype={'names': ('label', 'abund', 'abund2'),
                                                    'formats': (np.int, np.float, np.float)}).tolist()).transpose()

x = np.array([0.01, 0.1, 1, 10, 100, 1000])
y = np.array([a0001[2][len(a0001[2]) - 1], a01[2][len(a01[2]) - 1], a1[2][len(a1[2]) - 1], a10[2][len(a10[2]) - 1],
              a100[2][len(a100[2]) - 1], a1000[2][len(a1000[2]) - 1]])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.xlabel(r'Factor')
plt.ylabel(r"Mass fraction of {}^{25}Mg")
# print x, y
plt.plot(x, y, 'o-')
plt.xscale('log')

plt.show()
