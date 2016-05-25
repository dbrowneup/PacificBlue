#Written by Dan Browne

from pyfaidx import Fasta
from datetime import datetime
from string import maketrans
from math import ceil
import itertools as it
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


class ScaffoldGraph():

    def __init__(self, FastaFile, ConnectionLists, edge_cutoff=0.25):
        print "Entering ScaffoldGraph module:", str(datetime.now())
        self.G = nx.DiGraph()
        #Parse fasta sequences into graph as nodes
        print "1... Adding nodes to graph"
        seqs = Fasta(FastaFile)
        for x in seqs:
            self.G.add_node(int(x.name), seq=str(x))
            self.G.add_node(-1 * int(x.name), seq=self.revc(str(x)))
        print "Nodes in graph:", len(self.G.nodes())
        #Parse connections into graph as edges
        print "2... Adding edges to graph"
        for clist in ConnectionLists:
            for connect in clist:
                v1, v2, d, n = connect
                self.set_edge(v1, v2, d, n)
        print "Number of unfiltered edges:", self.G.size()
        self.degree_counter()
        #Filter noisy edges
        print "3... Removing noisy edges"
        for v1, v2 in self.G.edges():
            attr = zip(self.G[v1][v2]['W'], self.G[v1][v2]['D'])
            self.filter_noisy_edges(v1, v2, attr, edge_cutoff)
        print "Number of edges remaining:", self.G.size()
        self.degree_counter()
        #Filter weak edges from noisy nodes
        print "4... Removing weak edges"
        weak_edges = self.find_weak_edges()
        self.G.remove_edges_from(weak_edges)
        print "Number of edges remaining:", self.G.size()
        self.degree_counter()
        #Filter self-loop edges
        print "5... Removing self-loop edges"
        selfloops = self.find_selfloops()
        self.G.remove_edges_from(selfloops)
        print "Number of edges remaining:", self.G.size()
        self.degree_counter()
        #Test edge consistency and add edge complements
        print "6... Testing edge consistency and adding edge complements"
        comp_edges_added = 0
        for v1, v2, d in self.G.edges_iter(data=True):
            if self.G.has_edge(-1 * v2, -1 * v1):
                self.test_edge_consistency(v1, v2)
            else:
                self.G.add_edge(-1 * v2, -1 * v1, d)
                comp_edges_added += 1
        print "Number of edge complements added:", comp_edges_added
        print "Number of edges remaining:", self.G.size() 
        #Filter repetetive edges
        print "7... Filtering repetitive edges"
        repetitive_edges = self.find_repetitive_edges()
        self.G.remove_edges_from(repetitive_edges)
        print "Number of edges remaining:", self.G.size()
        self.degree_counter()
        #Filter edges leading to and from tips
        print "8... Removing tip edges"
        tip_edges = self.find_tip_edges()
        self.G.remove_edges_from(tip_edges)
        print "Number of edges remaining:", self.G.size()
        self.degree_counter()
        #Filter transitive edges from graph
        print "9... Removing transitive edges"
        transitive_edges = self.find_transitive_edges()
        self.G.remove_edges_from(transitive_edges)
        print "Number of edges remaining:", self.G.size()
        self.degree_counter()
        #Filter edges from branching nodes
        print "10... Removing branching edges"
        branching_edges = self.find_branches()
        self.G.remove_edges_from(branching_edges)
        print "Number of edges remaining:", self.G.size()
        self.degree_counter()
        print "Leaving ScaffoldGraph module:", str(datetime.now())
    
    def revc(self, seq):
        tr = maketrans('ATGCatgc', 'TACGtacg')
        return seq.translate(tr)[::-1]

    def set_edge(self, v1, v2, d, n):
        #W = weights: inversely proportional to number of connections
        #             reported from each PacBio read (n)
        #D = distance estimates: calculated from PacBio alignments
        try:
            self.G[v1][v2]['W'].append(1.0 / n)
            self.G[v1][v2]['D'].append(d)
        except:
            self.G.add_edge(v1, v2, {'W': [1.0 / n], 'D': [d]})

    def degree_counter(self):
        deg = self.G.in_degree().values()
        in_deg = [deg.count(n) for n in range(5)]
        plus = len(deg) - sum(in_deg)
        in_deg.append(plus)
        in_deg.append(max(deg))
        deg = self.G.out_degree().values()
        out_deg = [deg.count(n) for n in range(5)]
        plus = len(deg) - sum(out_deg)
        out_deg.append(plus)
        out_deg.append(max(deg))
        print "Degree Frequency:"
        print "IN  0:{0} 1:{1} 2:{2} 3:{3} 4:{4} 5+:{5} Max:{6}".format(*in_deg)
        print "OUT 0:{0} 1:{1} 2:{2} 3:{3} 4:{4} 5+:{5} Max:{6}".format(*out_deg)
    
    def filter_noisy_edges(self, v1, v2, attr, edge_cutoff):
        credible_attr = []
        for a in attr:
            if a[0] >= edge_cutoff:
                credible_attr.append(a)
        if len(credible_attr) == 0:
            self.G.remove_edge(v1, v2)
            return
        wt, de = zip(*credible_attr)
        if sum(wt) < 1:
            self.G.remove_edge(v1, v2)
            return
        de = np.array(de)
        if (max(de) - min(de)) > abs(de.mean()):
            self.G.remove_edge(v1, v2)
        else:
            self.G[v1][v2]['W'] = int(ceil(sum(wt)))
            self.G[v1][v2]['D'] = int(ceil(de.mean()))
    
    def find_weak_edges(self):
        weak_edges = []
        for n in self.G:
            if self.G.in_degree(n) > 2:
                in_wt = np.array([self.G[v1][v2]['W'] for v1, v2 in self.G.in_edges(n)])
                for v1, v2 in self.G.in_edges(n):
                    e_wt = self.G[v1][v2]['W']
                    if e_wt < in_wt.mean():
                        weak_edges.append((v1, v2))
            if self.G.out_degree(n) > 2:
                out_wt = np.array([self.G[v1][v2]['W'] for v1, v2 in self.G.out_edges(n)])
                for v1, v2 in self.G.out_edges(n):
                    e_wt = self.G[v1][v2]['W']
                    if e_wt < out_wt.mean():
                        weak_edges.append((v1, v2))
        return weak_edges
    
    def find_selfloops(self):
        selfloops = set(self.G.selfloop_edges())
        for e in self.G.edges_iter():
            v1, v2 = e
            if abs(v1) == abs(v2) and e not in selfloops:
                selfloops.add(e)
        return list(selfloops)
    
    def test_edge_consistency(self, v1, v2):
        de_f = self.G[v1][v2]['D']
        de_r = self.G[-1*v2][-1*v1]['D']
        de = np.array([de_f, de_r])
        if (max(de) - min(de)) > abs(de.mean()):
            self.G.remove_edge(v1, v2)
            self.G.remove_edge(-1*v2, -1*v1)
        else:
            self.G[v1][v2]['D'] = int(ceil(de.mean()))
            self.G[-1*v2][-1*v1]['D'] = int(ceil(de.mean()))
    
    def find_repetitive_edges(self):
        repetitive_edges = set([])
        for n in self.G:
            p = self.G.predecessors(n)
            s = self.G.successors(n)
            if len(p) < 2 or len(s) < 2:
                continue
            possibilities = it.product(p, s)
            for u, w in possibilities:
                if u == w:
                    continue
                if self.G.has_edge(u, w):
                    repetitive_edges.add((u, n))
                    repetitive_edges.add((n, w))
        return list(repetitive_edges)            
    
    def find_tip_edges(self):
        tip_edges = []
        for v1, v2 in self.G.edges_iter():
            v1_in = self.G.in_degree(v1)
            v1_out = self.G.out_degree(v1)
            v2_in = self.G.in_degree(v2)
            v2_out = self.G.out_degree(v2)
            if v1_out > 1 and v2_in == 1 and v2_out == 0:
                tip_edges.append((v1, v2))
            elif v1_out == 1 and v1_in == 0 and v2_in > 1:
                tip_edges.append((v1, v2))
        return tip_edges
    
    def find_transitive_edges(self):
        transitive_edges = []
        for n in self.G:
            p = self.G.predecessors(n)
            s = self.G.successors(n)
            if len(p) > 1 or len(p) == 0 or len(s) > 1 or len(s) == 0:
                continue
            if self.G.out_degree(p[0]) == 2 and self.G.in_degree(s[0]) == 2:
                if self.G.has_edge(p[0], s[0]):
                    transitive_edges.append((p[0], s[0]))
        return transitive_edges

    def find_branches(self):
        branching_edges = []
        for n in self.G:
            if len(self.G.in_edges(n)) > 1:
                branching_edges += self.G.in_edges(n)
            if len(self.G.out_edges(n)) > 1:
                branching_edges += self.G.out_edges(n)
        return branching_edges

