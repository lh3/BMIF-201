#!/usr/bin/env python

import random, math

# e[t][x]=e(x|t), emission probability of generating x from state t
e0 = [[0.5, 0.5], [0.9, 0.1]]

# p[s][t]=p(t|s) transition probability of s->t
p = [[0.9, 0.1], [0.2, 0.8]]

# q[t]=q(t) stationary distribution by solving linear equation P*x=x
q = [p[1][0]/(p[0][1]+p[1][0]), p[0][1]/(p[0][1]+p[1][0])]

def hmm_gen(l, p, e): # l: length; p: transition matrix; e: emission matrix
	t, y, x = 0, '', ''
	for i in range(l):
		r = random.random()
		if r < e[t][0]: a = 0
		elif r < e[t][0] + e[t][1]: a = 1
		else: a = 2
		y += "FC"[t]
		x += str(a)
		r = random.random()
		t = 0 if r < p[t][0] else 1
	return y, x

def hmm_forward(x, q, p, e): # x: observation; q: stationary; p: transition; e: emission
	L = len(x)
	f = [[0.0, 0.0] for i in range(L)]
	s = [1.0] * L
	a = int(x[0])
	f[0][0] = e[0][a] * q[0]
	f[0][1] = e[1][a] * q[1]
	for i in range(1, L):
		a = int(x[i])
		f[i][0] = e[0][a] * (f[i-1][0] * p[0][0] + f[i-1][1] * p[1][0])
		f[i][1] = e[1][a] * (f[i-1][0] * p[0][1] + f[i-1][1] * p[1][1])
		s[i] = f[i][0] + f[i][1]
		f[i][0] /= s[i]
		f[i][1] /= s[i]
	return f, s
	
def hmm_backward(x, q, p, e, s):
	L = len(x)
	b = [[0.0, 0.0] for i in range(L)]
	b[L-1] = [1, 1]
	for i in range(L - 2, -1, -1):
		a = int(x[i+1])
		b[i][0] = (b[i+1][0] * e[0][a] * p[0][0] + b[i+1][1] * e[1][a] * p[0][1]) / s[i+1]
		b[i][1] = (b[i+1][0] * e[0][a] * p[1][0] + b[i+1][1] * e[1][a] * p[1][1]) / s[i+1]
	return b

def hmm_em1(x, q, p, e00):
	e = [[0.5, 0.5], [e00, 1 - e00]]
	f, s = hmm_forward(x, q, p, e)
	b = hmm_backward(x, q, p, e, s)
	L = len(x)
	n0, n1, logP = 0, 0, 0
	for i in range(L):
		z = [f[i][0]*b[i][0], f[i][1] * b[i][1]]
		z0 = z[1] / (z[0] + z[1]) # actually z[0]+z[1] should always be 1
		logP += math.log(s[i])
		if int(x[i]) == 0:
			n0 += z0
		if int(x[i]) == 1:
			n1 += z0
	print("{:.2f}\t{:.2f}\t{:.4f} -> {:.4f}\t{}".format(n0, n1, e00, n0/(n0+n1), logP))
	return n0/(n0+n1)

# main function
random.seed(1)
y, x = hmm_gen(10000, p, e0)
e00 = 0.5
for i in range(50):
	e00 = hmm_em1(x, q, p, e00)
