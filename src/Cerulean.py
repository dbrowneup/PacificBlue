# This software is Copyright 2013 The Regents of the University of
# California. All Rights Reserved.
#
# Permission to copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without fee,
# and without a written agreement is hereby granted, provided that the above
# copyright notice, this paragraph and the following three paragraphs appear
# in all copies.
#
# Permission to make commercial use of this software may be obtained by
# contacting:
# Technology Transfer Office
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# invent@ucsd.edu
#
# This software program and documentation are copyrighted by The Regents of the
# University of California. The software program and documentation are supplied
# "as is", without any accompanying services from The Regents. The Regents does
# not warrant that the operation of the program will be uninterrupted or
# error-free. The end-user understands that the program was developed for
# research purposes and is advised not to rely exclusively on the program for
# any reason.
#
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO
# ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING
# OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
# THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
# CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
# ENHANCEMENTS, OR MODIFICATIONS.

#Author: Viraj Deshpande
#Contact: vdeshpan@eng.ucsd.edu

from long_contig_graph import LongContigGraph
from illumina_graph import Graph
from PacbioMapping import PacbioMapping
from PacbioAlignment import PacbioAlignment
from utils import Nx
from path.path import path
from sets import Set
from math import ceil
from datetime import datetime
import argparse

parser = argparse.\
ArgumentParser(description="Cerulean hybrid assembler v0.1.1")

parser.add_argument('--dataname', dest='dataName',
                    help="name of dataset", metavar='String',
                    action='store', type=str, nargs=1)

parser.add_argument('--basedir', dest='dataDir',
                    help="path to directory containing data", metavar='Path',
                    action='store', type=str, nargs=1)

parser.add_argument('--nproc', dest='nproc', default=[1],
                    help='number of parallel threads[1]', metavar='Integer',
                    action='store', type=int, nargs=1)

parser.add_argument('--mvl', dest='mvl', default=[2000],
                    help='minimum vertex length (default = 2000)', metavar='Integer',
                    action='store', type=int, nargs=1)
print 'Starting Cerulean.py:', str(datetime.now())
args = parser.parse_args()
dataName = args.dataName[0]
dataDir = args.dataDir[0]
nproc = args.nproc[0]
min_vertex_length = args.mvl[0]
prefix = dataDir + '/' + dataName
contigFa = prefix + '-scaffolds.fa'
contigDot = prefix + '-scaffolds.dot'
pbMapping = prefix + '_pacbio_contigs_mapping.fasta.m4'
pbMappingFormat = 'm4'
outFa = prefix + '-cerulean.fasta'
outDot = prefix + '-cerulean.dot'

print 'Beginning to load graph:', str(datetime.now())
ig = Graph()
ig.load(contigDot, contigFa)
print 'Graph loaded:', str(datetime.now())
print 'Number of vertices =', len(ig.vs)
print 'Number of edges =', len(ig.es)
original_contig_lengths = [v.length for v in ig.vs.values() if v.vid > 0]
original_contig_lengths.sort()
print ("Total input contig length ="), sum(original_contig_lengths)
original_contig_lengths.sort()
print "Input N50 =", Nx(original_contig_lengths)

print 'Beginning to load mapping:', str(datetime.now())
pbm = PacbioMapping(pbMapping, pbMappingFormat)
print "Mapping loaded:", str(datetime.now())
print "Number of mapping reads =", len(pbm.readToContig)
print "Number of mapping contigs =", len(pbm.contigToRead)
readlens = [PacbioAlignment(pbm.readToContig[r][0], pbMappingFormat).queryLen for r in pbm.readToContig]
readlens.sort()
print "Length of mapping reads =", sum(readlens)
print "N50 of mapping reads =", Nx(readlens)
pacbioCoverage = sum(readlens) / sum(original_contig_lengths)
print "Mapping Pacbio Coverage =", pacbioCoverage
min_read_threshold = ceil(pacbioCoverage * 3 / 17)
max_read_threshold = ceil(pacbioCoverage * 10 / 17)
print "Read Thresholds =", min_read_threshold, max_read_threshold

print 'Beginning to load LongContigGraph:', str(datetime.now())
lcg = LongContigGraph(ig, pbm, num_threads=nproc,
                      read_thresholds=(min_read_threshold,
                                       max_read_threshold))
print "LongContigGraph loaded:", str(datetime.now())
print 'Vertices:', len(lcg.vs) 
print 'Edges:', len(lcg.es)

print 'Beginning to process simple paths and splcg:', str(datetime.now())
splcg = lcg.get_simple_paths()
svs = Set([])
svs2 = Set([])
for p in splcg:
    for v in [p[0].v1, p[-1].v2]:
        if v.Name() in svs:
            svs2.add(v.Name())
    for e in p:
        svs.add(e.v1.Name())
        svs.add(e.v2.Name())
print 'Done processing simple paths and splcg:', str(datetime.now())

f = path(outDot).open('w')
f.write(lcg.__repr__())
f.close()

f = path(outFa).open('w')

seen_vertices = {vid: False for vid in ig.vs.keys()}

p_count = max(ig.vs.keys()) + 1
extended_contig_lengths = []
for p in splcg:
    s_start = 0
    s_end = len(p) - 1

    pathseq = ig.path_seq(p[s_start: s_end + 1], fill_gaps=False)
    if len(pathseq) < min_vertex_length:
        continue
    pathrep = ''
    for e in p[s_start: s_end + 1]:
        pathrep += str(e.v1.vid) + ';' + str(e.ovl) + ';'
    pathrep += str(e.v2.vid)

    f.write('>' + dataName + '|ref')
    fid = ['0'] * 7
    pcc = p_count
    for i in range(7):
        fid[i] = str(pcc % 10)
        pcc = pcc / 10
    f.write(''.join(fid[::-1]) + '\n')
    f.write(pathseq + '\n')
    extended_contig_lengths.append(len(pathseq))

    for e in p[s_start: s_end + 1]:
        seen_vertices[e.v1.vid] = True
        seen_vertices[e.v1.conj.vid] = True
        seen_vertices[e.v2.vid] = True
        seen_vertices[e.v2.conj.vid] = True
    p_count += 1

allc = splcg[:]

small_contig_length = 0
for vid in ig.vs.keys():
    if seen_vertices[vid]:
        continue
    if vid < 0:
        continue
    if ig.vs[vid].length < min_vertex_length:
        small_contig_length += ig.vs[vid].length
        continue
    f.write('>' + dataName + '|ref')
    fid = ['0'] * 7
    pcc = vid
    for i in range(7):
        fid[i] = str(pcc % 10)
        pcc = pcc / 10
    f.write(''.join(fid[::-1]) + '\n')
    f.write(pathseq + '\n')
    f.write(ig.vs[vid].seq + '\n')
    extended_contig_lengths.append(len(ig.vs[vid].seq))
    allc.append(ig.vs[vid])

print("# output contigs = "), len(extended_contig_lengths)
print ("Total output contig length = "), sum(extended_contig_lengths)
extended_contig_lengths.sort()
print "Output N50 = ", Nx(extended_contig_lengths)

print "Length of ignored contigs", small_contig_length
f.close()
