#!/usr/bin/env python

import random

# e[t][x]=e(x|t), emission probability of generating x from state t
e = [[0.5, 0.5], [0.9, 0.1]]

# p[s][t]=p(t|s) transition probability of s->t
p = [[0.9, 0.1], [0.2, 0.8]]

# q[t]=q(t) stationary distribution by solving linear equation P*x=x
q = [p[1][0]/(p[0][1]+p[1][0]), p[0][1]/(p[0][1]+p[1][0])]

def hmm_gen(l):
  t, ss, sl = 0, '', ''
  for i in range(l):
    r = random.random()
    if r < e[t][0]: x = 0
    elif r < e[t][0] + e[t][1]: x = 1
    else: x = 2
    ss += "FC"[t]
    sl += str(x)
    r = random.random()
    t = 0 if r < p[t][0] else 1
  return ss, sl

def hmm_forward(sl):
  L = len(sl)
  f = [[0.0, 0.0] for i in range(L)]
  x = int(sl[0])
  f[0][0] = e[0][x] * q[0]
  f[0][1] = e[1][x] * q[1]
  for i in range(1, L):
    x = int(sl[i])
    f[i][0] = e[0][x] * (f[i-1][0] * p[0][0] + f[i-1][1] * p[1][0])
    f[i][1] = e[1][x] * (f[i-1][0] * p[0][1] + f[i-1][1] * p[1][1])
  return f
  
def hmm_backward(sl):
  L = len(sl)
  b = [[0.0, 0.0] for i in range(L)]
  b[L-1][0] = 1
  b[L-1][1] = 1
  for i in range(L - 2, -1, -1):
    x = int(sl[i+1])
    b[i][0] = b[i+1][0] * e[0][x] * p[0][0] + b[i+1][1] * e[1][x] * p[0][1]
    b[i][1] = b[i+1][0] * e[0][x] * p[1][0] + b[i+1][1] * e[1][x] * p[1][1]
  return b

def hmm_decode(f, b):
  L, sd = len(f), ''
  for i in range(L):
    x = [f[i][0]*b[i][0], f[i][1] * b[i][1]]
    x0 = x[0] / (x[0] + x[1])
    l = 0 if x0 > 0.5 else 1
    print("{}\t{}\t{}".format(i, x0, x[0] + x[1]))
    sd += "FC"[l]
  return sd

random.seed(1)
ss, sl = hmm_gen(100)
print("Label:      {}".format(sl))
print("Truth:      {}".format(ss))

f = hmm_forward(sl)
b = hmm_backward(sl)
sd = hmm_decode(f, b)
print("Prediction: {}".format(sd))
