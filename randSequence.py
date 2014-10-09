

from collections import deque
import random

random.seed(0)

myseq = [random.choice(['a', 'c', 'g', 't']) for i in range(1000)]

print myseq[:10]       #this is 10

readlen = 25
starts = [random.randint(0, len(myseq)-readlen) for i in range(500)]
reads = [myseq[start:start+readlen] for start in  starts]

print reads[0]



results = "coverage.png"
coverage = [0] * len(myseq)
for start in starts:
  for i in range(start, start+readlen):
    coverage[i] += 1

import matplotlib, numpy
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(4,2))
plt.plot(coverage)
fig.tight_layout()
plt.savefig(results)
results


