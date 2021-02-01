#Our cohomological conditions gave rise to certain piecewise linear functions f on certain compact regions (themselves defined by finitely many linear conditions). The code below verifies the minimum values of each f. Since the functions are piecewise linear a lower bound is provided by the values of f at all points where a collection of boundary conditions intersect in a single point. Total run time ~ 1 min 30 sec

import itertools

#TETRAGONAL (Lemma 5.6, complement of B^\circ): verifies min(f) = 1/4 in specified region

def f((x1,x2,x3,y1,y2)):
	return 2*x3 - 2*x1 + y2 - y1

#inequalities in the definition of our region: ([a_1, a_2, a_3, b_1, b_2], c) means ax_1 + a_2x_2 + a_3x_3 + b_1y_1 + b_2y_2 \leq c
#the condition that x_1 + x_2 + x_3 = y_1 + y_2 = 1 will be imposed in next step
D1 = (([-1, 0, 0, 0, 0],0), ([1, -1, 0, 0, 0],0), ([0, 1, -1, 0, 0],0), ([0, 0, 0, 1, -1],0), ([2, 0, 0, 0, -1],0))
D2 = ()

def is_in_region(v): #returns True if vector v satisfies the inequalities specifed by D1
	for A in D1:
		if vector(A[0])*v > A[1]:
			return False
	return True

for A in itertools.combinations(D1+D2,3): #select three boundary conditions to impose (in addition to x_1 + x_2 + x_3 = y_1 + y_2 = 1)
	M = Matrix([[1, 1, 1, 0, 0], [0, 0, 0, 1, 1], A[0][0], A[1][0], A[2][0]])
	if M.det() != 0: #because f is piecewise linear, it suffices to consider a collection of independent boundary conditions
		v = vector([1, 1, A[0][1], A[1][1], A[2][1]])
		sol = M.inverse()*v #the point determined by the intersection of these boundary conditions
		if is_in_region(sol):
			if f(sol) < 1/4:
				print 'f is less than 1/4 at ', sol

print f((1/4,3/8,3/8,1/2,1/2)) == 1/4

#TETRAGONAL (Lemma 5.8, complement of H^circ): verifies min(f) = 1/4 in specified region

def f((x1,x2,x3,y1,y2)):
	return 2*x3 - 2*x1 + y2 - y1 - max(0, y2 - 2*x1) - max(0, y2 - x1 - x2)

#inequalities in the definition of our region: ([a_1, a_2, a_3, b_1, b_2], c) means ax_1 + a_2x_2 + a_3x_3 + b_1y_1 + b_2y_2 \leq c
#the condition that x_1 + x_2 + x_3 = y_1 + y_2 = 1 will be imposed in next step
D1 = (([-1, 0, 0, 0, 0],0), ([1, -1, 0, 0, 0],0), ([0, 1, -1, 0, 0],0), ([0, 0, 0, -1, 0],0), ([-2, 0, 0, 1, 0],0), ([2, 0, 0, 0, -1],0), ([0, -2, 0, 0, 1],0), ([-1, 0, -1, 0, 1],0))

#additional boundary conditions where the function f changes
D2 = (([-2, 0, 0, 0, 1], 0), ([-1, -1, 0, 0, 1], 0))

for A in itertools.combinations(D1+D2,3): #select three boundary conditions to impose (in addition to sum of coords = 1)
	M = Matrix([[1, 1, 1, 0, 0], [0, 0, 0, 1, 1], A[0][0], A[1][0], A[2][0]])
	if M.det() != 0:
		v = vector([1, 1, A[0][1], A[1][1], A[2][1]])
		sol = M.inverse()*v #the point determined by the intersection of these boundary conditions
		if is_in_region(sol):
			if f(sol) < 1/4:
				print 'f is less than 1/4 at ', sol

print f((1/4,3/8,3/8,1/2,1/2)) == 1/4

#PENTAGONAL (Lemma 5.12, complement of B^\circ): verifies min(f) = 1/5 in specified region

def f((x1,x2,x3,x4,y1,y2,y3,y4,y5)):
	return 3*x4 + x3 - x2 - 3*x1 + 4*y5 + 2*y4 - 2*y2 - 4*y1

#inequalities in the definition of our region: ([a_1, a_2, a_3, a_4, b_1, b_2, b_3, b_4, b_5], c) means ax_1 + ... + a_4x_4 + b_1y_1 + ... + b_5y_5 \leq c
#the condition that x_1 + ... + x_4 = 1 and y_1 + ... + y_5 = 2 will be imposed in the next step
D1 = (([-1, 0, 0, 0, 0, 0, 0, 0, 0], 0), ([1, -1, 0, 0, 0, 0, 0, 0, 0], 0), ([0, 1, -1, 0, 0, 0, 0, 0, 0], 0), ([0, 0, 1, -1, 0, 0, 0, 0, 0], 0), ([0, 0, 0, 0, -1, 0, 0, 0, 0], 0), ([0, 0, 0, 0, 1, -1, 0, 0, 0], 0), ([0, 0, 0, 0, 0, 1, -1, 0, 0], 0), ([0, 0, 0, 0, 0, 0, 1, -1, 0], 0), ([0, 0, 0, 0, 0, 0, 0, 1, -1], 0), ([1, 0, 0, 0, 1, 1, 0, 0, 0], 1))

D2 = ()

for A in itertools.combinations(D1+D2,7): #select seven boundary conditions to impose (in addition to x_1 + ... + x_4 = 1 and y_1 + ... + y_5 = 2)
	M = Matrix([[1, 1, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1, 1], A[0][0], A[1][0], A[2][0], A[3][0], A[4][0], A[5][0], A[6][0]])
	if M.det() != 0:
		v = vector([1, 2, A[0][1], A[1][1], A[2][1], A[3][1], A[4][1], A[5][1], A[6][1]])
		sol = M.inverse()*v #the point determined by the intersection of these boundary conditions
		if is_in_region(sol):
			if f(sol) < 1/5:
				print 'f is less than 1/5 at ', sol

print f((1/5,4/15,4/15,4/15,2/5,2/5,2/5,2/5,2/5)) == 1/5


#PENTAGONAL (Lemma 5.13, complement of H^circ): verifies min(f) = 1/5 in specified region

def f((x1,x2,x3,x4,y1,y2,y3,y4,y5)):
	return 3*x4 + x3 - x2 - 3*x1 + 4*y5 + 2*y4 - 2*y2 - 4*y1 - sum([max(0, 1 - y1 - y2 - x) for x in [x1, x2, x3, x4]]) - sum([max(0, 1 - y1 - y3 - x) for x in [x1, x2, x3]]) - sum([max(0, 1 - y1 - y4 - x) for x in [x1, x2]]) - sum([max(0, 1 - y2 - y3 - x) for x in [x1, x2]])

#inequalities in the definition of our region: ([a_1, a_2, a_3, a_4, b_1, b_2, b_3, b_4, b_5], c) means ax_1 + ... + a_4x_4 + b_1y_1 + ... + b_5y_5 \leq c
#the condition that x_1 + ... + x_4 = 1 and y_1 + ... + y_5 = 2 will be imposed in the next step
D1 = (([-1, 0, 0, 0, 0, 0, 0, 0, 0], 0), ([1, -1, 0, 0, 0, 0, 0, 0, 0], 0), ([0, 1, -1, 0, 0, 0, 0, 0, 0], 0), ([0, 0, 1, -1, 0, 0, 0, 0, 0], 0), ([0, 0, 0, 0, -1, 0, 0, 0, 0], 0), ([0, 0, 0, 0, 1, -1, 0, 0, 0], 0), ([0, 0, 0, 0, 0, 1, -1, 0, 0], 0), ([0, 0, 0, 0, 0, 0, 1, -1, 0], 0), ([0, 0, 0, 0, 0, 0, 0, 1, -1], 0), ([0, 0, 0, -1, -1, 0, -1, 0, 0], -1), ([0, 0, -1, 0, -1, 0, 0, -1, 0], -1), ([0, 0, -1, 0, 0, -1, -1, 0, 0], -1), ([1, 0, 0, 0, 1, 1, 0, 0, 0], 1))

#additional boundary conditions where the function f changes:
D2 = (([1, 0, 0, 0, 1, 1, 0, 0, 0], 1), ([0, 1, 0, 0, 1, 1, 0, 0, 0], 1), ([0, 0, 1, 0, 1, 1, 0, 0, 0], 1), ([0, 0, 0, 1, 1, 1, 0, 0, 0], 1), ([1, 0, 0, 0, 1, 0, 1, 0, 0], 1), ([0, 1, 0, 0, 1, 0, 1, 0, 0], 1), ([0, 0, 1, 0, 1, 0, 1, 0, 0], 1), ([1, 0, 0, 0, 1, 0, 0, 1, 0], 1), ([0, 1, 0, 0, 1, 0, 0, 1, 0], 1), ([1, 0, 0, 0, 0, 1, 1, 0, 0], 1), ([0, 1, 0, 0, 0, 1, 1, 0, 0], 1))

for A in itertools.combinations(D1+D2,7): #select seven boundary conditions to impose (in addition to  x_1 + ... + x_4 = 1 and y_1 + ... + y_5 = 2)
	M = Matrix([[1, 1, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1, 1], A[0][0], A[1][0], A[2][0], A[3][0], A[4][0], A[5][0], A[6][0]])
	if M.det() != 0:
		v = vector([1, 2, A[0][1], A[1][1], A[2][1], A[3][1], A[4][1], A[5][1], A[6][1]])
		sol = M.inverse()*v #the point determined by the intersection of these boundary conditions
		if is_in_region(sol):
			if f(sol) < 1/5:
				print 'f is less than 1/5 at ', sol

print f((1/5,4/15,4/15,4/15,2/5,2/5,2/5,2/5,2/5)) == 1/5


