#Programmers: H.Ortiz-Zuazaga, O.Rosado Ramirez, R.Lopez Torres
#
#
#This program is meant to solve the sequence assebmbly problem by
#mixing all the known possible solutions. It is called CopperHead.
#


import networkx as nx
import random
#import matplotlib, numpy
import find_euler_path
from graphviz import Digraph
import matplotlib.pyplot as plt




G = nx.Graph()

DG=nx.DiGraph()

random.seed(0)
myseq = "".join([random.choice(['a', 'c', 'g', 't']) for i in range(1000)]) #RANGE ES 1000
print myseq[:60], "..."

readlen = 25
starts = [random.randint(0, len(myseq)-readlen) for i in range(500)]  #range es 500
reads = [myseq[start:start+readlen] for start in  starts]
print reads[0]

def makebad(read):
    point = random.randint(0,len(read)-1)
    bad = read[:point] + "x" + read[point+1:]
    return bad

badreads = 15
for i in range(badreads):
    reads[-i] = makebad(reads[-i])

print reads    #[-20:]  #humberto hizo este print no se pa que

print '\n' # a~ade un espacio ahi para hacerme sentir mejor

#print reads #este animalito es el que vamos a despedazar ahora

#-------------------------------Function by: O.Rosado 2015------------------------


def SequenceDNA(): #TURN THIS INTO FUNCTION AND USE MAP
    print "----------------------------------"
    print "Esto es lo nuevo"

    k = 3

    kmers = []

    KmerList = []


#-------------------------function by R.Lopez Torres----------------------------
#This function creates the kmers from the strands.


    for read in reads:

        for i in range(len(read) - k):

            kmer = read[i:i + k]

            kmers.insert(i, kmer)


        KmerList.insert(i, kmers)

        kmers = []


#-------------------------function by R.Lopez Torres----------------------------

    print KmerList

    print "\n"

    print "The length of the kmer list is: ", len(KmerList)

    print "\n"

#this function removes the kmers that least appear
    for kmer in KmerList:

        for k in kmer:

            if kmer.count(k) < 5:

                kmer.remove(k)
#this function removes the kmers that least appear

#Aqui empieza la funcion que crea el grafo q no esta haciendo na!!
    for kmers in KmerList:
        #for k in kmers:

            #G.add_node(k) si le pongo este, necesita el segundo for
        G.add_path(kmers)

    #G.nx.Graph(m)
    #nx.draw(G)    #esta estaba funcionando pero hizo el grafo gey ese
    #nx.draw_networkx(G)
    #nx.Graph()
    #draw_nx(m)
    nx.write_dot(G,"kmers.dot")  #esto no hace na Humberto!
    #p = nx.spring_layout( G )

    #nx.draw(G,p)

    #plt.savefig( "g.png" )

    print G.number_of_nodes()
    print G.nodes()


"""
    for k in Kmers:

        Kmers.count(k)

        #print "\n"

    #for kr in range(len(Kmers)):

        if Kmers.count(k) < 5:
            Kmers.remove(k)


    print Kmers

    print "\n"

    return Kmers
"""

SequenceDNA()


#find_euler_path.find_eulerian_path

#THIS IS UNCALLED. MUST EDIT LATER
def AssembleDNA(Kmers):

    newSeq = Kmers[0]

    for i in range(1, len(Kmers)):
        newSeq += Kmers[i][-1]


    #for k in range(len(Kmers)):
    #    ker = ''.join(Kmers[0] + Kmers[k:])

    print "this is the assembly\n"



    print newSeq

    print "\n"

    print myseq


    if newSeq in myseq:
        print "NOBEL PRIZE"
    else:
        print "The strands arrent the same. :,("


    return newSeq

#THIS IS UNCALLED. MUST FIX LATER

#AssembleDNA(Kmers)


#this Sequences the sequences into kmers

"""

m = 0

for m in range(len(reads)):
    if reads[i] == reads[i + 1]:
        del list[i]
    else:
        m += 1

print "number of reads is ", m
"""
#esto es lo raro que esto y haciendo

"""
ESTO VERIFICA SI LOS KMERS TIENEN READS FALSOS.
err = 'x'
for i in range(len(reads)):
    if err in reads[i]:
        print "heh"
    else:
        print "hoh"
"""

print "es hasta aqui"
print "-----------------------------------"

#def termino(newSeq):


#termino(newSeq)








#------------------------------------end of Function by: O.Rosado--------------

"""
coverage = [0] * len(myseq)
for start in starts:
  for i in range(start, start+readlen):
    coverage[i] += 1


fig=plt.figure(figsize=(4,2))
plt.plot(coverage)
fig.tight_layout()

k = 9

for read in reads:
  for i in range(len(read) - k - 1):
    kmer_a = read[i:i+k]
    kmer_b = read[i+1:i+k+1]
    DG.add_edge(kmer_a, kmer_b)

components = list(nx.connected_components(DG.to_undirected()))
len(components)

map(len, components)
#components.close()

placement = [0] * len(myseq)

for i in range(len(components)):
    for read in components[i]:
        try:
            place = myseq.index(read)
        except ValueError:
            continue
        placement[place] = i

placefig=plt.figure(figsize=(4,2))
plt.plot(placement)
placefig.tight_layout()

best = nx.subgraph(DG, components[0])

def eulerian_path(G, source=None):
    #Find an Eulerian path in a weakly connected component, stolen from eulerian_circuit
    if not nx.is_weakly_connected(G):
        raise nx.NetworkXError("G is not connected.")

    g = G.__class__(G) # copy graph structure (not attributes)

    # set starting node
    if source is None:
        for v in g.nodes():
            indeg = g.in_degree(v)
            outdeg = g.out_degree(v)
            if outdeg == (indeg + 1):
                break
        else:
            raise nx.NetworkXError("G doesn't have an Eulerian circuit.")
    else:
        v = source

    while g.out_degree(v) > 0:
        n = v
        bridge = None
        # sort nbrs here to provide stable ordering of alternate cycles
        nbrs = sorted([v for u,v in g.edges(n)])
        for v in nbrs:
            g.remove_edge(n,v)
            bridge = not nx.is_connected(g.to_undirected())
            if bridge:
                g.add_edge(n,v)  # add this edge back and try another
            else:
                break  # this edge is good, break the for loop
        if bridge:
            g.remove_edge(n,v)
            g.remove_node(n)
        yield (n,v)


path = eulerian_path(best)
nodepath = [u for u,v in path]
seq = nodepath[0] + "".join([node[-1] for node in nodepath[1:]])
seq in myseq
len(seq)

#-------------------------Function by: O.Rosado-----2015-------------------------------------
def oldSkin():
    print "==========================================================="
    print "This is the sequence with errors. Errors are represented with X's."
    print seq
    print "===========================================================\n"


oldSkin()
"""#NUEVO HUM


"""   COMENTANDO SHEDDING
Sequence = []  #Making an empty list, home of the kstring

#Shedding function makes remover artificial errors from the string
def Shedding():
    print "----------------------------------------------------"
    print "This is the sequence without errors."
    for j in range(len(seq)):
        str = seq[j]
        Sequence.append (str)
        if 'x' in Sequence:
            Sequence.remove('x')

    print ''.join(Sequence)
    print "----------------------------------------------------\n"


Shedding()

#--------------------------Function by O.Rosado----------2015=---------------------------
"""










#"---------------------------NO NECESITO ESTO--------------------------------"


#parte de la grafica que no voy a usar ahora
"""
i = 0
for contig in components:
    cg = nx.subgraph(DG, contig)
    nx.write_dot(cg, "contig-" + str(i) + "-with-errors.dot")
    i += 1

nx.write_dot(DG, "full-graph.dot")

from IPython.display import SVG
SVG(filename='contig-0.svg')
"""
