#Written by Dan Browne

from pyfaidx import Fasta
from datetime import datetime
from string import maketrans
from math import ceil
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


class ScaffoldGraph():

    def __init__(self, FastaFile, ConnectionLists, edge_cutoff=0.25):
        print "Entering ScaffoldGraph module:", str(datetime.now())
        self.graph = nx.DiGraph()
        #Parse fasta sequences into graph as nodes
        print "Adding nodes to graph"
        seqs = Fasta(FastaFile)
        for x in seqs:
            self.graph.add_node(int(x.name), seq=str(x))
            self.graph.add_node(-1 * int(x.name), seq=self.revc(str(x)))
        print "Nodes in graph:", len(self.graph.nodes())
        #Parse connections into graph as edges
        print "Adding edges to graph"
        for clist in ConnectionLists:
            for connect in clist:
                v1, v2, d, n = connect
                self.set_edge(v1, v2, d, n)
        print "Number of unfiltered edges:", len(self.graph.edges())
        self.degree_counter()
        #Filter noisy edges
        print "Removing noisy edges"
        for v1, v2 in self.graph.edges():
            attr = zip(self.graph[v1][v2]['weights'], self.graph[v1][v2]['dist_estimates'])
            self.test_edge(v1, v2, attr, edge_cutoff)
        print "Number of edges remaining:", len(self.graph.edges())
        self.degree_counter()
        #Filter self-loop edges
        print "Removing self-loop edges"
        selfloops = self.graph.selfloop_edges()
        self.graph.remove_edges_from(selfloops)
        print "Number of edges remaining:", len(self.graph.edges())
        self.degree_counter()
        #Filter edges leading to tips
        print "Removing tip edges"
        tip_edges = self.find_tips()
        self.graph.remove_edges_from(tip_edges)
        print "Number of edges remaining:", len(self.graph.edges())
        self.degree_counter()
        #Filter weak edges from noisy nodes
        print "Removing weak edges"
        weak_edges = self.find_weak_edges()
        self.graph.remove_edges_from(weak_edges)
        print "Number of edges remaining:", len(self.graph.edges())
        self.degree_counter()
        #Filter edges from branching nodes
        print "Removing branching edges"
        branching_edges = self.find_branches()
        self.graph.remove_edges_from(branching_edges)
        print "Number of edges remaining:", len(self.graph.edges())
        self.degree_counter()
        print "Leaving ScaffoldGraph module:", str(datetime.now())
    
    def revc(self, seq):
        tr = maketrans('ATGCatgc', 'TACGtacg')
        return seq.translate(tr)[::-1]

    def set_edge(self, v1, v2, d, n):
        try:
            self.graph[v1][v2]['weights'].append(float(1) / n)
            self.graph[v1][v2]['dist_estimates'].append(d)
        except:
            self.graph.add_edge(v1, v2, {'weights': [float(1) / n], 'dist_estimates': [d]})

    def degree_counter(self):
        deg = self.graph.in_degree().values()
        in_deg = [deg.count(n) for n in range(5)]
        plus = len(deg) - sum(in_deg)
        in_deg.append(plus)
        in_deg.append(max(deg))
        deg = self.graph.out_degree().values()
        out_deg = [deg.count(n) for n in range(5)]
        plus = len(deg) - sum(out_deg)
        out_deg.append(plus)
        out_deg.append(max(deg))
        print "Degree Frequency:"
        print "IN 0:{0} 1:{1} 2:{2} 3:{3} 4:{4} 5+:{5} Max:{6}".format(*in_deg)
        print "OUT 0:{0} 1:{1} 2:{2} 3:{3} 4:{4} 5+:{5} Max:{6}".format(*out_deg)
    
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
            self.graph[v1][v2]['weights'] = int(sum(wt))
            self.graph[v1][v2]['dist_estimates'] = int(ceil(de.mean()))
    
    def find_weak_edges(self):
        weak_edges = []
        for n in self.graph.nodes():
            in_wt = np.array([self.graph[v1][v2]['weights'] for v1, v2 in self.graph.in_edges(n)])
            out_wt = np.array([self.graph[v1][v2]['weights'] for v1, v2 in self.graph.out_edges(n)])
            for v1, v2 in self.graph.in_edges(n):
                e_wt = self.graph[v1][v2]['weights']
                if e_wt < in_wt.mean():
                    weak_edges.append((v1, v2))
            for v1, v2 in self.graph.out_edges(n):
                e_wt = self.graph[v1][v2]['weights']
                if e_wt < out_wt.mean():
                    weak_edges.append((v1, v2))
        return weak_edges
    
    def find_tips(self):
        tip_edges = []
        for v1, v2 in self.graph.edges():
            v1_out = self.graph.out_degree(v1)
            v2_in = self.graph.in_degree(v2)
            v2_out = self.graph.out_degree(v2)
            if v1_out > 1 and v2_in == 1 and v2_out == 0:
                tip_edges.append((v1, v2))
        return tip_edges

    def find_branches(self):
        branching_edges = []
        for n in self.graph.nodes():
            if len(self.graph.in_edges(n)) > 1:
                branching_edges += self.graph.in_edges(n)
            if len(self.graph.out_edges(n)) > 1:
                branching_edges += self.graph.out_edges(n)
        return branching_edges

