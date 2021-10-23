def em1(q, l, a): # one round of EM
  N, R = len(q), len(a) # #transcripts, #reads
  if len(q) != len(l): raise "Error!"
  r = [0 for t in range(N)] # this keeps the next round of q[]
  for i in range(R):
    s = 0.0
    for t in a[i]:
      if t < N:
        s += q[t] / l[t]  # compute the sum
    if s == 0.0: continue # avoid division by zero
    for t in a[i]:
      if t < N:
        r[t] += q[t] / l[t] / s
  for t in range(N):
    q[t] = r[t] / R;

def q2p(q, l): # convert nucleotide fraction q to transcript fraction p
  N = len(q)
  p = [q[t] / l[t] for t in range(N)]
  s = sum(p)
  p = [p[t] / s for t in range(N)]
  return p

l = [200, 400, 150]
a = [[0,1], [0,1], [1], [1], [2]]

q = [1./len(l)] * len(l)
for i in range(25):
  em1(q, l, a)
  p = q2p(q, l)
  print("\t".join(["{0:.4f}".format(j) for j in p]))
