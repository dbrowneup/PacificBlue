#Re-written by Dan Browne on 04/29/16

from PacbioSubgraph import PacbioSubgraph
from pyfaidx import Fasta
from datetime import datetime
from string import maketrans
from multiprocessing import Pool
from math import ceil
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#A million thanks for this hack, found here: http://www.rueckstiess.net/research/snippets/show/ca1d7d90
def UnwrapParallelSubgraph(arg, **kwarg):
    return LongContigGraph.ParallelSubgraph(*arg, **kwarg)

class LongContigGraph():

    def __init__(self, FastaFile, PacbioMapping, num_threads=4, edge_cutoff=(float(1) / 4)):
        print "Beginning LongContigGraph:", str(datetime.now())
        self.graph = nx.DiGraph()
        self.mapping = PacbioMapping
        #Parse fasta sequences into graph as nodes
        seqs = Fasta(FastaFile)
        for x in seqs:
            self.graph.add_node(int(x.name), seq=str(x))
            self.graph.add_node(-1 * int(x.name), seq=self.revc(str(x)))
        #Process alignments in parallel
        sg_pool = Pool(processes=num_threads)
        reads = self.mapping.readToContig.keys()
        chunk = int(ceil(float(len(reads))/num_threads))
        connection_lists = sg_pool.map(func=UnwrapParallelSubgraph, iterable=zip([self]*len(reads), reads), chunksize=chunk)
        connection_lists = [clist for clist in connection_lists if clist is not None]
        for clist in connection_lists:
            for connect in clist:
                v1, v2, d, n = connect
                self.set_edge(v1, v2, d, n)
        #Filter noisy edges
        for v1, v2 in self.graph.edges():
            attr = zip(self.graph[v1][v2]['weights'], self.graph[v1][v2]['dist_estimates'])
            self.test_edge(v1, v2, attr, edge_cutoff)
        #Filter weak edges from noisy nodes
        for n in self.graph:
            self.weak_edges(n)
        #Filter edges from branching nodes
        for n in self.graph:
            self.filter_branches(n)
        print "Finishing LongContigGraph:", str(datetime.now())
    
    def ParallelSubgraph(self, read):
        sg = PacbioSubgraph(read, self.mapping)
        self.mapping.readToContig[read] = sg.mapping.readToContig[read]
        if len(sg.Connects) == 0:
            return
        reported_connects = []
        for a1, a2, d in sg.Connects:
            v1 = -1 * a1.tName if a1.tStrand == 1 else a1.tName
            v2 = -1 * a2.tName if a2.tStrand == 1 else a2.tName
            n = len(sg.Connects)
            reported_connects.append((v1, v2, d, n))
        return reported_connects
    
    def revc(self, seq):
        tr = maketrans('ATGCatgc', 'TACGtacg')
        return seq.translate(tr)[::-1]

    def set_edge(self, v1, v2, d, n):
        try:
            self.graph[v1][v2]['weights'].append(float(1) / n)
            self.graph[v1][v2]['dist_estimates'].append(d)
        except:
            self.graph.add_edge(v1, v2, {'weights': [float(1) / n], 'dist_estimates': [d]})

    def test_edge(self, v1, v2, attr, edge_cutoff):
        credible_attr = []
        for a in attr:
            if a[0] >= edge_cutoff:
                credible_attr.append(a)
        if len(credible_attr) == 0:
            self.graph.remove_edge(v1, v2)
            return
        wt, de = zip(*credible_attr)
        if sum(wt) < 1:
            self.graph.remove_edge(v1, v2)
            return
        de = np.array(de)
        if (max(de) - min(de)) > abs(de.mean()):
            self.graph.remove_edge(v1, v2)
        else:
            self.graph[v1][v2]['weights'] = sum(wt)
            self.graph[v1][v2]['dist_estimates'] = de.mean()
    
    def weak_edges(self, n):
        in_wt = np.array([self.graph[v1][v2]['weights'] for v1, v2 in self.graph.in_edges(n)])
        out_wt = np.array([self.graph[v1][v2]['weights'] for v1, v2 in self.graph.out_edges(n)])
        for v1, v2 in self.graph.in_edges(n):
            e_wt = self.graph[v1][v2]['weights']
            if e_wt < in_wt.mean():
                self.graph.remove_edge(v1, v2)
        for v1, v2 in self.graph.out_edges(n):
            e_wt = self.graph[v1][v2]['weights']
            if e_wt < out_wt.mean():
                self.graph.remove_edge(v1, v2)
    
    def filter_branches(self, n):
        if len(self.graph.in_edges(n)) > 1:
            self.graph.remove_edges_from(self.graph.in_edges(n))
        if len(self.graph.out_edges(n)) > 1:
            self.graph.remove_edges_from(self.graph.out_edges(n))

