#Written by Dan Browne on 05/11/16

from datetime import datetime
import networkx as nx
import numpy as np
import matplotlib
import maplotlib.pyplot as plt
import sys

class PathFinder():
    
    def __init__(self, scaffold_graph):
        self.G = scaffold_graph
        #Build strandless list of sequences
        sequences = set([x for x in self.G.nodes() if x > 0])
        #Define weakly connected components
        component_graphs = set([x for x in nx.weakly_connected_component_subgraphs(self.G)])
        single_node_graphs = set([x for x in component_graphs if len(x.nodes()) == 1])
        multi_node_graphs = set([x for x in component_graphs if len(x.nodes()) > 1])
        #Define paths in multi-node graphs
        seqs_in_multi_node_paths = set([])
        multi_node_paths = [nx.dag_longest_path(x) for x in multi_node_graphs]
        for path in multi_node_paths:
            seqs_in_multi_node_paths |= set(path)
        #Remove single node seqs complementary to multi-node seqs
        for g in single_node_graphs:
            if -1 * g.nodes()[0] in seqs_in_multi_node_paths:
                single_node_graphs.discard(g)
        #Build multi-node paths into scaffolds
        self.scaffolds = []
        for g in multi_node_graphs:
            self.build_scaffold(g)
        #Add single node seqs to scaffolds list

    def build_scaffold(self, graph):
        scaffold_order = []
        scaffold_seq = ''
        #Add sequence and edge data to scaffold_order
        for v, e in zip(graph.nodes(), graph.edges() + [None]):
            seq = graph[v]['seq']
            try:
                dist = graph[v][e[1]]['dist_estimates']
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
                scaffold_seq = self.merge_sequences(merge_seq, scaffold_seq, item)
            else:
                sys.exit('FATAL ERROR: Unknown scaffold building operation!')
        #Report assembled scaffold sequence
        self.scaffolds.append(scaffold_seq)

    def merge_sequences(self, merge_seq, scaff_seq, dist_est):
        five_prime = np.array([n for n in merge_seq])
        thre_prime = np.array([n for n in scaff_seq])
        #Test for matching bases in overlap region
        boolean_overlap = five_prime[dist_est:] == thre_prime[:abs(dist_est)]
        #Report matching bases, mask mismatches
        sequence_overlap = (five_prime[dist_est + i] if boolean_overlap[i] else 'N' for i in range(abs(dist_est)))
        sequence_overlap = ''.join(sequence_overlap)
        #Return full_sequence
        full_sequence = ''.join(five_prime[:dist_est]) + sequence_overlap + ''.join(thre_prime[abs(dist_est):])
        return full_sequence

