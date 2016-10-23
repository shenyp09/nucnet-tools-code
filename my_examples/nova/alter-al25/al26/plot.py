from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


p1 = np.array(np.loadtxt('p1.txt', dtype={'names': ('label', 'time'),
                                          'formats': (np.int, np.float)}).tolist()).transpose()

p10 = np.array(np.loadtxt('p10.txt', dtype={'names': ('label', 'time'),
                                            'formats': (np.int, np.float)}).tolist()).transpose()

p100 = np.array(np.loadtxt('p100.txt', dtype={'names': ('label', 'time'),
                                              'formats': (np.int, np.float)}).tolist()).transpose()

p1000 = np.array(np.loadtxt('p1000.txt', dtype={'names': ('label', 'time'),
                                                'formats': (np.int, np.float)}).tolist()).transpose()

a0001 = np.array(np.loadtxt('al26_0.001.dat', dtype={'names': ('label', 'abund'),
                                                     'formats': (np.int, np.float)}).tolist()).transpose()
a01 = np.array(np.loadtxt('al26_0.1.dat', dtype={'names': ('label', 'abund'),
                                                 'formats': (np.int, np.float)}).tolist()).transpose()
a1 = np.array(np.loadtxt('al26_1.dat', dtype={'names': ('label', 'abund'),
                                              'formats': (np.int, np.float)}).tolist()).transpose()
a10 = np.array(np.loadtxt('al26_10.dat', dtype={'names': ('label', 'abund'),
                                                'formats': (np.int, np.float)}).tolist()).transpose()
a100 = np.array(np.loadtxt('al26_100.dat', dtype={'names': ('label', 'abund'),
                                                  'formats': (np.int, np.float)}).tolist()).transpose()
a1000 = np.array(np.loadtxt('al26_1000.dat', dtype={'names': ('label', 'abund'),
                                                    'formats': (np.int, np.float)}).tolist()).transpose()

m1 = np.ones(len(a1[1]))
m10 = 10 * np.ones(len(a10[1]))
m100 = 100 * np.ones(len(a100[1]))
m1000 = 1000 * np.ones(len(a1000[1]))

x = np.array([0.01, 0.1, 1, 10, 100, 1000])
y = np.array([a0001[1][len(a0001[1]) - 1], a01[1][len(a01[1]) - 1], a1[1][len(a1[1]) - 1], a10[1][len(a10[1]) - 1],
              a100[1][len(a100[1]) - 1], a1000[1][len(a1000[1]) - 1]])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.xlabel(r'Factor')
plt.ylabel(r"Mass fraction of {}^{26}Al")
# print x, y
plt.plot(x, y, 'o-')
# plt.yscale('log')
plt.xscale('log')

plt.show()
