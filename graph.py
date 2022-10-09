import numpy as np
from numpy import cos, pi
import scipy as sp

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271]
P = len(primes)

### todo: make this be part of a proper AbelianGroup class
def groupSum(t1, t2, group):
	assert len(t1) == len(t2) and len(t1) == len(group), "groupSum length mismatch"
	n = len(t1)
	s = [ 0 for j in range(n) ]
	for j in range(n):
		s[j] = (t1[j] + t2[j]) % primes[group[j][0]]**group[j][1]
	return tuple(s)
assert groupSum((19,),(11,),[(3,1)]) == (2,), "groupSum error"

def lexValue(t, group):
	assert len(t) == len(group), "lexValue length mismatch"
	n = len(t)
	tot = 0
	for j in range(n):
		modulus = primes[group[j][0]]**group[j][1]
		assert 0 <= t[j] < modulus, "not a properly-formed group element"
		if j > 0:
			tot *= modulus
		tot += t[j]
	return tot

assert lexValue((0,0), [(3,1),(4,1)]) == 0, "lexValue error"
assert lexValue((4,5), [(3,1),(4,1)]) == 49, "lexValue error"

class AbelianCayley:
	
	def __init__(self, group, generators, aMat = True, aList = False, givenAMat = None):
		### todo: enforce that group generators are sorted first by prime then by exponent, enforce generators to be symmetric, enforce that generators are properly-formed group elements, make inputs be using an AbelianGroup class, enforce that generators has no duplicates
		self.rank = len(group)
		self.order = 1
		for (i,e) in group:
			assert i < P, "Prime too large"
			self.order *= primes[i]**e
		self.group = group
		self.group_elements = [ [n] for n in range(primes[self.group[0][0]]**self.group[0][1]) ]
		for j in range(1, self.rank):
			(i, e) = self.group[j]
			while len(self.group_elements[0]) <= j:
				t = self.group_elements.pop(0)
				for n in range(primes[i]**e):
					self.group_elements.append(t + [n])
		self.group_elements = [ tuple(t) for t in self.group_elements ] # self.group_elements comes endowed with a natural (lexicographic) order
		self.gens = generators
		self.deg = len(generators)
		
		self.hasAMat = False
		if aMat or not (givenAMat is None):
			self.hasAMat = True

		if not (givenAMat is None):
			self.adjMat = givenAMat
		elif aMat:
			self.makeAMat()
		
		self.find_spectrum()
		
		### todo: adjacency list
		
	def makeAMat(self):
		self.adjMat = np.zeros((self.order, self.order), dtype=int)
		for x in self.group_elements:
			# print(x)
			Lx = lexValue(x, self.group)
			for s in self.gens:
				xs = groupSum(x, s, self.group)
				Lxs = lexValue(xs, self.group)
				# print(x,Lx,s,xs,Lxs)
				self.adjMat[Lx, Lxs] = 1
		self.hasAMat = True
	
	def times(self, G2, wantAMat = True):
		r1, r2 = self.rank, G2.rank
		o1, o2 = self.order, G2.order
		g1, g2 = self.group, G2.group
		gen1, gen2 = self.gens, G2.gens
		if not wantAMat:
			return AbelianCayley( g1 + g2, [e + (0,)*r2 for e in gen1] + [(0,)*r1 + e for e in gen2], False, False, None)
		elif self.hasAMat and G2.hasAMat:
			a1, a2 = self.adjMat, G2.adjMat
			return AbelianCayley( g1 + g2, [e + (0,)*r2 for e in gen1] + [(0,)*r1 + e for e in gen2], False, False, np.kron(a1, np.eye(o2, dtype='int')) + np.kron(np.eye(o1, dtype='int'), a2) ) ### todo: order of kron may be wrong...
		else:
			AbelianCayley( g1 + g2, [e + (0,)*r2 for e in gen1] + [(0,)*r1 + e for e in gen2], True, False, None)
	
	def __mul__(self, G2):
		return self.times(G2)
	def __rmul__(self, G2):
		return self.times(G2)
	
		
	### to implement: str?, format?, repr?, eq?, ne?
	
	def find_spectrum(self):
		self.spectrum = []
		for x in self.group_elements:
			t = 0
			for s in self.gens:
				p = 0
				for i in range(self.rank):
					p += x[i]*s[i]/primes[self.group[i][0]]**self.group[i][1]
				t += cos(2*pi*p)
			self.spectrum.append(t/self.deg)
		self.spectrum = np.sort(1 - np.array(self.spectrum))
		self.gap = self.spectrum[1]
		self.connected = (self.gap > 0)
		

"""
GC = AbelianCayley([(3, 1)], [(1,), (6,)])
print(GC.order, GC.rank)
print(GC.group_elements)
print(GC.adjMat)
"""
"""
GC2 = AbelianCayley([(1, 1), (2, 1)], [(1, 0), (0, 1), (2, 0), (0, 4)])
print(GC2.order, GC2.rank)
print(GC2.group_elements)
print(GC2.adjMat)
"""

"""
# check adjMat for product graphs
GC2 = AbelianCayley([(0, 1)], [(1,)])
# print(GC2.adjMat)
GC3 = AbelianCayley([(1, 1)], [(1,), (2,)])
# print(GC3.adjMat)
GC23 = GC2*GC3
print(GC23.group, GC23.gens)
print(GC23.adjMat)
GCalt = AbelianCayley(GC23.group, GC23.gens)
print(GCalt.adjMat)
"""

"""
# check spectrum and gap
GC2 = AbelianCayley([(0, 1)], [(1,)])
print(GC2.spectrum, GC2.gap)
GC3 = AbelianCayley([(1, 1)], [(1,), (2,)])
print(GC3.spectrum, GC3.gap)
GC23 = GC2*GC3
print(GC23.spectrum, GC23.gap)
"""

"""
# check connectivity
GC4 = AbelianCayley([(0, 2)], [(2,)])
print(GC4.connected)
"""
