#Written by Dan Browne on 05/11/16

from datetime import datetime
from string import maketrans
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
        print "1... Defining weakly connected components"
        component_graphs = set([g for g in nx.weakly_connected_component_subgraphs(self.G)])
        single_node_graphs = set([g for g in component_graphs if len(g.nodes()) == 1])
        multi_node_graphs = set([g for g in component_graphs if len(g.nodes()) > 1])
        print "Number of single-node components:", len(single_node_graphs)
        print "Number of multi-node components:", len(multi_node_graphs)
        #Consolidate unscaffolded nodes, discard reverse strand
        print "2... Consolidating single-node components"
        unscaffolded = set([g.nodes()[0] for g in single_node_graphs])
        discard_nodes = set([n for n in unscaffolded if n < 0])
        for g in iter(single_node_graphs.copy()):
            if g.nodes()[0] in discard_nodes:
                single_node_graphs.discard(g)
        print "Number of unscaffolded sequences:", len(single_node_graphs)
        #Classify multi-node graphs
        print "3... Classifying multi-node components"
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
        #Build scaffolds from DAGs
        print "4... Building scaffolds from directed acyclic graphs"
        self.scaffolds = set([])
        for g in DAG:
            self.build_dag_scaffold(g)
        #Consolidating complementary scaffolds, keep first found
        print "5... Consolidating complementary scaffolds"
        consolidated_scaff = set([])
        for seq in iter(self.scaffolds):
            comp = self.revc(seq)
            if comp in self.scaffolds:
                if comp not in consolidated_scaff:
                    consolidated_scaff.add(seq)
            else:
                print "WARNING: non-complemented scaffold"
        self.scaffolds = consolidated_scaff
        print "Number of scaffolds assembled:", len(self.scaffolds)
        #Build scaffolds from Eulerian graphs
        
        #Add unscaffolded seqs to scaffolds list
        print "6... Adding unscaffolded sequences to output"
        for g in single_node_graphs:
            seq = self.G.node[g.nodes()[0]]['seq']
            self.scaffolds.add(seq)
        print "Leaving PathFinder module:", str(datetime.now())
    
    def build_dag_scaffold(self, g):
        scaffold_order = []
        scaffold_seq = ''
        #Determine path in DAG
        path = nx.dag_longest_path(g)
        #Add sequence and edge data to scaffold_order
        for n in path:
            seq = g.node[n]['seq']
            try:
                s = next(g.successors_iter(n))
                dist = g[n][s]['D']
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
        self.scaffolds.add(scaffold_seq)
    
    def merge_seqs(self, merge_seq, scaff_seq, dist):
        #Return sequences joined with masked overlap
        five_prime = merge_seq[:dist]
        thre_prime = scaff_seq[abs(dist):]
        overlap = 'N' * abs(dist)
        full_sequence = five_prime + overlap + thre_prime
        return full_sequence
    
    def revc(self, seq):
        tr = maketrans('ATGCatgc', 'TACGtacg')
        return seq.translate(tr)[::-1]
