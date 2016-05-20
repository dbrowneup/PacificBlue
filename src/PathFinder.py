#Written by Dan Browne on 05/11/16

from datetime import datetime
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys

class PathFinder():
    
    def __init__(self, scaffold_graph):
        print "Entering PathFinder module:", str(datetime.now())
        self.G = scaffold_graph.copy()
        #Build strandless list of sequences
        sequences = set([n for n in self.G.nodes() if n > 0])
        #Define weakly connected components
        component_graphs = set([g for g in nx.weakly_connected_component_subgraphs(self.G)])
        single_node_graphs = set([g for g in component_graphs if len(g.nodes()) == 1])
        multi_node_graphs = set([g for g in component_graphs if len(g.nodes()) > 1])
        print "Number of single-node graphs:", len(single_node_graphs)
        print "Number of multi-node graphs:", len(multi_node_graphs)
        #Classify multi-node graphs
        DAG = set([])
        Euler = set([])
        for g in multi_node_graphs:
            if nx.is_directed_acyclic_graph(g):
                DAG.add(g)
            elif nx.is_eulerian(g):
                Euler.add(g)
            else:
                sys.exit("FATAL ERROR: Unknown multi-node graph type!")
        print "Number of directed acyclic graphs:",  len(DAG)
        print "Number of Eulerian graphs:", len(Euler)
        #Determine nodes in single-node and multi-node graphs
        self.multi_nodes = set([])
        for g in multi_node_graphs:
            self.multi_nodes |= set(g.nodes())
        self.single_nodes = set([g.nodes()[0] for g in single_node_graphs])
        #Consolidate all nodes
        self.unscaffolded = set([])
        self.single_scaffolded = set([])
        self.double_scaffolded = set([])
        for n in sequences:
            self.classify_nodes(n)
        single_node_graphs = [g for g in single_node_graphs if self.filter_sng(g)]
        #Build scaffolds from DAGs
        self.scaffolds = []
        for g in DAG:
            self.build_dag_scaffold(g)
        #Build scaffolds from Eulerian graphs
        
        #Add single node seqs to scaffolds list
        for g in single_node_graphs:
            seq = self.G.node[g.nodes()[0]]['seq']
            self.scaffolds.append(seq)
        print "Leaving PathFinder module:", str(datetime.now())
    
    def classify_nodes(self, n):
        if n in self.single_nodes and -1*n in self.single_nodes:
            self.unscaffolded.add(n)
        elif n in self.single_nodes and -1*n in self.multi_nodes:
            self.single_scaffolded.add(-1*n)
        elif n in self.multi_nodes and -1*n in self.single_nodes:
            self.single_scaffolded.add(n)
        elif n in self.multi_nodes and -1*n in self.multi_nodes:
            self.double_scaffolded.add(n)
            self.double_scaffolded.add(-1*n)
    
    def filter_sng(self, g):
        n = g.nodes()[0]
        if n in self.unscaffolded:
            return True
        else:
            return False
    
    def build_dag_scaffold(self, g):
        scaffold_order = []
        scaffold_seq = ''
        #Determine path in DAG
        path = nx.dag_longest_path(g)
        #Add sequence and edge data to scaffold_order
        for n in path:
            seq = self.G.node[n]['seq']
            try:
                e = g.out_edges(n)[0]
                dist = g[n][e[1]]['D']
                scaffold_order.append(seq)
                scaffold_order.append(dist)
            except:
                scaffold_order.append(seq)
        #Convert gaps in scaffold_order into strings of Ns
        for i in range(len(scaffold_order)):
            value = scaffold_order[i]
            if type(value) is int and value > 0:
                scaffold_order[i] = 'N' * value
        #Build scaffold_seq and merge overlaps
        while len(scaffold_order) > 0:
            merge_seq = ''
            item = scaffold_order.pop()
            if type(item) is str:
                scaffold_seq = item + scaffold_seq
            elif type(item) is int and item == 0:
                merge_seq = scaffold_order.pop()
                scaffold_seq = merge_seq + scaffold_seq
            elif type(item) is int and item < 0:
                merge_seq = scaffold_order.pop()
                scaffold_seq = self.merge_seqs(merge_seq, scaffold_seq, item)
            else:
                sys.exit('FATAL ERROR: Unknown scaffold building operation!')
        #Report assembled scaffold sequence
        self.scaffolds.append(scaffold_seq)

    def merge_seqs(self, merge_seq, scaff_seq, dist):
        #Return sequences joined with masked overlap
        five_prime = merge_seq[:dist]
        thre_prime = scaff_seq[abs(dist):]
        overlap = 'N' * abs(dist)
        full_sequence = five_prime + overlap + thre_prime
        return full_sequence

