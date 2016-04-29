#Re-written by Dan Browne on 04/29/16

from PacbioSubgraph import PacbioSubgraph
from utils import parmap, rc
from pyfaidx import Fasta
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math


class LongContigGraph():

    def __init__(self, FastaFile, PacbioMapping, num_threads=1):
        self.graph = nx.Graph()
        self.mapping = PacbioMapping
        self.num_threads = num_threads
        mapping_reads = self.mapping.readToContig.keys()
        batch_size = int(math.ceil(float(len(mapping_reads)) / self.num_threads)) + 1
        print "Number of threads = ", self.num_threads
        print "Batch size: ", batch_size
        #Parse fasta sequences into graph as nodes
        seqs = Fasta(FastaFile)
        for x in seqs:
            self.graph.add_node(int(x.name), seq=str(x))
            self.graph.add_node(-1 * int(x.name), seq=rc(str(x)))
        #Process alignments in parallel
        def parallel_subgraphs(thread_id=0):
            slice_start = thread_id * batch_size
            slice_end = min((thread_id + 1) * batch_size, len(mapping_reads))
            for r in mapping_reads[slice_start:slice_end]:
                sg = PacbioSubgraph(r, self.mapping)
                if len(sg.Connects) == 0:
                    continue
                for x in sg.Connects:
                    v1 = -1 * x[0].tName if x[0].tStrand == 1 else x[0].tName
                    v2 = -1 * x[1].tName if x[1].tStrand == 1 else x[1].tName
                    dist = x[2]
                    print "v1:", v1, "v2:", v2, "dist:", dist, "weight", 1/len(sg.Connects)
                    try:
                        self.graph[v1][v2]['weight'] += (float(1) / len(sg.Connects))
                        self.graph[v1][v2]['dist_estimates'].append(dist)
                        print "Edge already existed!"
                        return
                    except:
                        self.graph.add_edge(v1, v2, {'weight': (float(1) / len(sg.Connects)), 'dist_estimates': [dist]})
                        print "Created new edge"
                        return
        parmap(parallel_subgraphs, range(self.num_threads))
