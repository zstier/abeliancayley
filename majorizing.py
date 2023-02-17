from graph import *
from random import randrange
import numpy as np
from numpy import sqrt

def majorize(a, b, prec=0): # what is the minimum index k such that for j < k, a[j] <= b[j], to additive error eps?
	# inputs must be sorted, increasing
	for k in range(min(len(a),len(b))):
		if a[k] + prec > b[k]:
			return k
	return min(len(a),len(b))
assert majorize([0, 0.25, 0.5],[0, 0.3, 0.3]) == 2, "majorize error"


# group data
(p,e) = (3,1) # (3,1) corresponds to Z/7Z
base = primes[p]**e
d = 4 # power
order = base**d
print("order of group is", order)
stdgens = []
for i in range(d):
	v, vinv = [0 for _ in range(d)], [0 for _ in range(d)]
	v[i] = 1
	vinv[i] = base-1
	stdgens.append(tuple(v))
	stdgens.append(tuple(vinv))
# print(stdgens)

G = AbelianCayley([(p,e) for _ in range(d)], stdgens)
# print(G.spectrum)

T = 1000 # number of trials
data = []
for _ in range(T):
	a = 3 # number of additional pairs of generators. could be made random, if desired
	newgens = []
	for x in stdgens:
		newgens.append(x)
	for g in range(a):
		v = [randrange(base) for _ in range(d)]
		tv = tuple(v)
		if tv in newgens or tv == (0,)*d:
			g -= 1 # accidentally repeated a generator (either a standard one or a random one already generated), so re-roll
		else:
			vinv = [(base-v[i])%base for i in range(d)] # v's inverse
			tvinv = tuple(vinv)
			newgens.append(tv)
			newgens.append(tvinv)
	H = AbelianCayley([(p,e) for _ in range(d)], newgens)
	# print(newgens)
	# print(H.spectrum)
	k = majorize(G.spectrum, H.spectrum)
	data.append(k)
	if k < sqrt(order):
		print(newgens), k
print(data)
print(min(data))
