#Re-written by Dan Browne on 04/29/16

from PacbioSubgraph import PacbioSubgraph
from pyfaidx import Fasta
from datetime import datetime
from string import maketrans
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


class LongContigGraph():

    def __init__(self, FastaFile, PacbioMapping, edge_cutoff=1):
        print "Beginning LongContigGraph:", str(datetime.now())
        self.graph = nx.DiGraph()
        self.mapping = PacbioMapping
        #Parse fasta sequences into graph as nodes
        seqs = Fasta(FastaFile)
        for x in seqs:
            self.graph.add_node(int(x.name), seq=str(x))
            self.graph.add_node(-1 * int(x.name), seq=self.revc(str(x)))
        #Process alignments
        for r in self.mapping.readToContig.keys():
            sg = PacbioSubgraph(r, self.mapping)
            self.mapping.readToContig[r] = sg.mapping.readToContig[r]
            if len(sg.Connects) == 0:
                continue
            for x in sg.Connects:
                v1 = -1 * x[0].tName if x[0].tStrand == 1 else x[0].tName
                v2 = -1 * x[1].tName if x[1].tStrand == 1 else x[1].tName
                dist = x[2]
                try:
                    self.graph[v1][v2]['weights'].append(float(1) / len(sg.Connects))
                    self.graph[v1][v2]['dist_estimates'].append(dist)
#                    print "Edge already existed!"
                except:
                    self.graph.add_edge(v1, v2, {'weights': [float(1) / len(sg.Connects)], 'dist_estimates': [dist]})
#                    print "Created new edge"
        #Filter noisy edges
        for e in self.graph.edges():
            v1, v2 = e
            wt, de = zip(self.graph[v1][v2]['weights'], self.graph[v1][v2]['dist_estimates'])
            self.test_edge(v1, v2, wt, de, edge_cutoff)
        print "Finishing LongContigGraph:", str(datetime.now())
    
    def revc(self, seq):
        tr = maketrans('ATGC', 'TACG')
        return seq.translate(tr)[::-1]

