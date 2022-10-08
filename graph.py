import numpy as np
import scipy as sp

""" # maybe useful later?
class AbelianGroup:
	
	def __init__(self,
"""

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
assert (groupSum((19,),(11,),[(3,1)]) == (2,)), "groupSum error"

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

assert (lexValue((0,0), [(3,1),(4,1)]) == 0), "lexValue error"
assert (lexValue((4,5), [(3,1),(4,1)]) == 49), "lexValue error"

class AbelianCayley:
	
	def __init__(self, group, generators, aMat = True, aList = False, givenAMat = None):
		### todo: enforce that group generators are sorted first by prime then by exponent, enforce generators to be symmetric, enforce that generators are properly-formed group elements, make inputs be using an AbelianGroup class, enforce that generators has no duplicates
		self.rank = len(group)
		self.order = 1
		for (i,e) in group:
			assert (i < P), "Prime too large"
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
		
		if givenAMat != None:
			self.adjMat = givenAMat
		elif aMat:
			self.adjMat = np.zeros((self.order, self.order), dtype=int)
			for x in self.group_elements:
				print(x)
				Lx = lexValue(x, self.group)
				for s in self.gens:
					xs = groupSum(x, s, self.group)
					Lxs = lexValue(xs, self.group)
					print(x,Lx,s,xs,Lxs)
					self.adjMat[Lx, Lxs] = 1
		
		### todo: adjacency list
		
		
	### to implement: str?, format?, repr?, eq?, ne?

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
