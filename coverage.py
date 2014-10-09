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