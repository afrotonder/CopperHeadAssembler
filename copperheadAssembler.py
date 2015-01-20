
#THIS IS ORG MODE
#Programmers: Humberto Ortiz-Zuazaga,  Omar Rosado Ramirez
#Description: This program is made to solve the main problems in sequence assembly.
#             The program, called CopperHead, is written in Python 2 and uses scientific libraries
#             such as matplotlib. This is the current working model, properly called CopperHead.2,
#             which takes on one of 2 existing (relevant to this current study) problems*.
#
# * Ortiz-Zuazaga, H., Rosado Ramirez, O.(2014). CopperHead: Rethinking Sequence Assembly.
#
#Comments: This program is in its infancy and does little of what is wanted. The rest will come
#          in the next version.
#
#

from collections import deque
import random

random.seed(0)

myseq = [random.choice(['a', 'c', 'g', 't']) for i in range(1000)]


print myseq[:10]

readlen = 25
starts = [random.randint(0, len(myseq)-readlen) for i in range(500)]
reads = [myseq[start:start+readlen] for start in  starts]

print reads[0]

#snaking begins here

for i in range(494,500):
  reads[i][random.randint(1, 10)] = 'x'        #Takes random x's in the given range and
print reads    #reads[-10:]                    #inserts them in the reads, simulating errors.

#snaking ends here

"""


#No tengo la libreria de matplotlib y no la necesito ahora.
#Esto esta comentado porque lo unico que hace es crear una grafica.
#Pregutar a Humberto despues si la grafica es data necesaria.

results = "coverage.png"    #This produces the graph.


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
"""

k = 16
kmers = []
kmergraph = {}

for read in reads:
  read = "".join(read)
  for i in range(len(read) - k - 1):
    kmer_a = read[i:i+k]
    kmer_b = read[i+1:i+k+1]
    if kmer_a not in kmers:
       kmers.append(kmer_a)
    if kmer_b not in kmers:
       kmers.append(kmer_b)
    if kmer_a not in kmergraph.keys():
       kmergraph[kmer_a] = deque([kmer_b])
    elif kmer_b not in kmergraph[kmer_a]:
       kmergraph[kmer_a].append(kmer_b)
    if kmer_b not in kmergraph.keys():
       kmergraph[kmer_b] = deque([])

print kmers[:10]


def findstart(kmers, kmergraph, path):
  start = None
  for i in range(len(path)):
    if len(kmergraph[kmers[path[i]]]) > 0:
      start = path[i]
      break
  return start

def extendpath(kmers, kmergraph, path):
  start = findstart(kmers, kmergraph, path)
  splicepoint = start
  if start != None:
    newpath = []
    while len(kmergraph[kmers[start]]) > 0:
      next_kmer = kmergraph[kmers[start]].popleft()
      start = kmers.index(next_kmer)
      newpath.append(start)
    path[splicepoint:splicepoint+1] = newpath
  return start


for i in range(330,331):
  start = i
  print start,
  path = [start]
  next = extendpath(kmers, kmergraph, path)
  while (next != None):
    next = extendpath(kmers, kmergraph, path)
  print len(path)
  k = 16
  kmers = []
  kmergraph = {}

  for read in reads:
    read = "".join(read)
    for i in range(len(read) - k - 1):
      kmer_a = read[i:i+k]
      kmer_b = read[i+1:i+k+1]
      if kmer_a not in kmers:
	 kmers.append(kmer_a)
      if kmer_b not in kmers:
	 kmers.append(kmer_b)
      if kmer_a not in kmergraph.keys():
	 kmergraph[kmer_a] = deque([kmer_b])
      elif kmer_b not in kmergraph[kmer_a]:
	 kmergraph[kmer_a].append(kmer_b)
      if kmer_b not in kmergraph.keys():
	 kmergraph[kmer_b] = deque([])

  print kmers[:10]


  assembly = kmers[path[0]]
for i in range(1, len(path)):
  assembly += kmers[path[i]][-1]

len(assembly)



def wrap(seq):
  start = 0
  for i in range(60, len(seq), 60):
    print seq[start: i]
    start = i
  print seq[start:]


wrap(assembly)

nx.write_dot(DG, "full-graph.dot") #Shows graph path with the Gephi program


#Estas son las notas, pensamientos e ideas que ocurren durante el experimento.
"""
oct/12/2014
I have two options. Make an if statement which snakes the sequence and removes the errors(x's).
What just occured to me is that i can make a similar if statement that removes the errors
from the reads instead of the sequences.

nov/2014
the above statement was made aprox two weeks before today, which is incorrect because the errors
must be removed in the reads, which is where they are created, or in snakng process, inserted.
"""
