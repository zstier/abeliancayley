from graph import AbelianCayley, P, primes
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt

# want to find graphs with large spectral gaps

gaps = []
for n in range(1, P):
	G = AbelianCayley([(n, 1)], [(1,), (primes[n]//2,)])
	G.find_spectrum()
	gaps.append(G.gap)

gaps = np.array(gaps)

"""
plt.plot([primes[i] for i in range(1, P)], gaps)
plt.plot([n for n in range(1, primes[P-1])], [1/sqrt(n) for n in range(1, primes[P-1])])
plt.show()
plt.close()
"""
plt.plot([primes[i] for i in range(1, P)], 1/sqrt(gaps))
plt.show()
# """
