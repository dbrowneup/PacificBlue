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

import sys
import argparse
from datetime import datetime
from multiprocessing import Pool
from math import ceil

import matplotlib
import matplotlib.pyplot as plt

from PacbioMapping import PacbioMapping
from PacbioSubgraph import PacbioSubgraph
from ScaffoldGraph import ScaffoldGraph
from PathFinder import PathFinder


#Parse command-line arguments

parser = argparse.ArgumentParser(description="PacificBlue Genome Scaffolder")

parser.add_argument('--blasr-alignments', help="BLASR alignment file", dest='blasr_alignments', metavar='FILE')
parser.add_argument('--blasr-format', help="Format of BLASR alignments", dest='blasr_format', default="m1", choices=['m1', 'm4', 'm5'])
parser.add_argument('--fasta', help="Fasta sequences to scaffold", dest='fasta_file', metavar='FILE')
parser.add_argument('--threads', help="Number of CPUs to use", dest='num_threads', metavar='NN', type=int, default=1)
parser.add_argument('--output', help="Name of output file", dest='output_file', metavar='FILE')
parser.add_argument('--cov_cutoff', help="PacbioSubgraph coverage cutoff (default=3)", dest='cov_cutoff', metavar='NN', type=int, default=3)
parser.add_argument('--fraction', help="PacbioSubgraph repeat fraction (default=0.5)", dest='fraction', metavar='(0,1]', type=float, default=0.5)
parser.add_argument('--edge_cutoff', help="ScaffoldGraph edge cutoff (default=0.25)", dest='edge_cutoff', metavar='(0,1]', type=float, default=0.25)

args = parser.parse_args()


#Worker function for parallel processing of read alignments

def ParallelSubgraph(item):
    read, mapping = item
    sg = PacbioSubgraph(read, mapping, cov_cutoff=args.cov_cutoff, fraction=args.fraction)
    if len(sg.Connects) == 0:
        return
    n = len(sg.Connects)
    report = []
    for a1, a2, d in sg.Connects:
        v1 = -1 * a1.tName if a1.tStrand == 1 else a1.tName
        v2 = -1 * a2.tName if a2.tStrand == 1 else a2.tName
        report.append((v1, v2, d, n))
    return report


#Main function to run PacificBlue pipeline

def main():
    #Load BLASR alignments into PacbioMapping object
    mapping = PacbioMapping(args.blasr_alignments, fileFormat=args.blasr_format)
    #Process read alignments in parallel
    print "Entering parallel PacbioSubgraph module:", str(datetime.now())
    print "Number of threads to use:", args.num_threads
    sg_pool = Pool(processes=args.num_threads)
    reads = mapping.readToContig.keys()
    chunk = int(ceil(float(len(reads))/args.num_threads))
    connection_lists = sg_pool.map(func=ParallelSubgraph, iterable=mapping.readToContig.items(), chunksize=chunk)
    sg_pool.close()
    sg_pool.join()
    connection_lists = [clist for clist in connection_lists if clist is not None]
    print "Leaving parallel PacbioSubgraph module:", str(datetime.now())
    #Load reported connections in ScaffoldGraph object
    scaff = ScaffoldGraph(args.fasta_file, connection_lists, edge_cutoff=args.edge_cutoff)
    #Build scaffold sequences with PathFinder object
    paths = PathFinder(scaff.G)
    #Write output to Fasta file
    out = open(args.output_file, "w")
    for i in range(len(paths.scaffolds)):
        out.write(">scaffold_"+str(i+1)+"\n")
        out.write(paths.scaffolds[i]+"\n")
    out.close()


if __name__== '__main__':
    main()
