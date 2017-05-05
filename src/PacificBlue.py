#!/usr/bin/env python

#PacificBlue Genome Scaffolding Tool for PacBio Long Reads
#
#Written by Dan Browne
#
#Devarenne Lab
#Department of Biochemistry & Biophysics
#Texas A&M University
#College Station, TX
#
#Contact: dbrowne.up@gmail.com

#Load modules

import argparse
import sys
from datetime import datetime
from math import ceil
from multiprocessing import Pool

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from PacbioMapping import PacbioMapping
from PacbioSubgraph import PacbioSubgraph
from ScaffoldGraph import ScaffoldGraph
from PathFinder import PathFinder


#Parse command-line arguments

parser = argparse.ArgumentParser(description="PacificBlue Genome Scaffolder")

parser.add_argument('--blasr-alignments', help="BLASR alignment file", \
                    dest='blasr_alignments', metavar='FILE')
parser.add_argument('--blasr-format', help="Format of BLASR alignments", \
                    dest='blasr_format', default="m1", choices=['m1', 'm4', 'm5'])
parser.add_argument('--fasta', help="Fasta sequences to scaffold", \
                    dest='fasta_file', metavar='FILE')
parser.add_argument('--threads', help="Number of CPUs to use", \
                    dest='num_threads', metavar='NN', type=int, default=1)
parser.add_argument('--chunksize', help="Number of tasks per process" \
                    dest='chunksize', metavar='INT', type=int, default=100)
parser.add_argument('--output', help="Name of output file", \
                    dest='output_file', metavar='FILE')
parser.add_argument('--cov_cutoff', help="PacbioSubgraph coverage cutoff (default=3)", \
                    dest='cov_cutoff', metavar='INT', type=int, default=3)
parser.add_argument('--fraction', help="PacbioSubgraph repeat fraction (default=0.5)", \
                    dest='fraction', metavar='(0,1]', type=float, default=0.5)
parser.add_argument('--edge_cutoff', help="ScaffoldGraph edge cutoff (default=0.25)", \
                    dest='edge_cutoff', metavar='(0,1]', type=float, default=0.25)
parser.add_argument('--edge_weight', help="ScaffoldGraph edge weight requirement (default=1)", \
                    dest='edge_weight', metavar='INT', type=int, default=1)

args = parser.parse_args()


#Worker function for parallel processing of read alignments

def subgraph(alignments):
    report = []
    with PacbioSubgraph(alignments, args.cov_cutoff, args.fraction) as sg:
        connects = sg.make_connects()
        if len(connects) == 0:
            return
        else:
            n = len(connects)
            for a1, a2, d in connects:
                v1 = -1 * a1.tName if a1.tStrand == 1 else a1.tName
                v2 = -1 * a2.tName if a2.tStrand == 1 else a2.tName
                report.append((v1, v2, d, n))
    return report


#Main function to run PacificBlue pipeline

def main():
    #Load BLASR alignments into PacbioMapping object
    print "Entering parallel PacbioSubgraph module:", str(datetime.now())
    print "Number of threads:", args.num_threads
    print "Size of data chunks:", args.chunk
    mapping = PacbioMapping(args.blasr_alignments, fileFormat=args.blasr_format, \
                            chunksize=args.chunksize, num_threads=args.num_threads)
    #Process read alignments in parallel
    connections = []
    counter = 0
    for chunk in mapping:
        sg_pool = Pool(processes=args.num_threads, maxtasksperchild=args.chunksize)
        for x in sg_pool.imap_unordered(subgraph, chunk.itervalues(), chunksize=args.chunksize):
            if x is not None:
                connections.append(x)
        sg_pool.close()
        sg_pool.join()
        counter += 1
    print "Finished processing "+str(counter)+" chunks"
    print "Leaving parallel PacbioSubgraph module:", str(datetime.now())
    #Load reported connections in ScaffoldGraph object
    scaff = ScaffoldGraph(args.fasta_file, connections, \
                        edge_cutoff=args.edge_cutoff, edge_weight=args.edge_weight)
    #Build scaffold sequences with PathFinder object
    paths = PathFinder(scaff.G)
    #Write output to Fasta file
    out = open(args.output_file, "w")
    for i, seq in enumerate(paths.scaffolds):
        out.write(">scaffold_"+str(i+1)+"\n")
        out.write(seq+"\n")
    out.close()


if __name__== '__main__':
    main()
