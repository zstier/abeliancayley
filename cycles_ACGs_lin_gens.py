from graph import AbelianCayley, P, primes
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, log2, floor
from sympy.ntheory import factorint


# want to find graphs with large spectral gaps

# """
gaps = []

max = 400

for n in range(5, max): # primes[P-1], 2):
	g = n//5
	d = factorint(n)
	ps = d.keys() # primes dividing n
	gp = [] # characterization of the group
	gens = [] # generators of the group. already includes first generator (1)
	new_gens = [[] for _ in range(1,g+1)] # corresp. to 1, ..., log n
	neg_new_gens = [[] for _ in range(1,g+1)] # corresp. to -1, ..., -log n
	for p in ps:
		gp.append((primes.index(p), d[p]))
		pe = p**d[p]
		for x in range(1,g+1):
			new_gens[x-1].append(x % pe)
			neg_new_gens[x-1].append((-x) % pe)
	gens += [tuple(x) for x in new_gens]
	gens += [tuple(x) for x in neg_new_gens]
	G = AbelianCayley(gp, gens)
	G.find_spectrum()
	gaps.append(G.gap)

gaps = np.array(gaps)

# plt.plot([n for n in range(3, 1259, 2)], gaps)
plt.plot([n for n in range(5, max)], gaps)
plt.show()
plt.close()
plt.plot([n for n in range(5, max)], 1/sqrt(gaps))
plt.show()
plt.close()

"""

gaps = []

for n in primes[1:]: # primes[P-1], 2):
	d = factorint(n)
	ps = d.keys() # primes dividing n
	gp = [] # characterization of the group
	gens = [(1,)*len(ps)] # generators of the group. already includes first generator (1)
	new_gen = [] # corresp. to n//2. will be converted (cast) into tuple
	for p in ps:
		gp.append((primes.index(p), d[p]))
		pe = p**d[p]
		new_gen.append( (n//2) % pe )
	gens.append(tuple(new_gen))
	# if n < 100: print(gp, gens)
	G = AbelianCayley(gp, gens)
	G.find_spectrum()
	gaps.append(G.gap)

gaps = np.array(gaps)

print(gaps)

"""
