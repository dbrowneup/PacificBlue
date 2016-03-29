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

#Author: Eric Fung
#Contact: kspham@eng.ucsd.edu

import experimental
import random
# from utils import plot
from loadgraph import loadgraph
from PacbioMapping import PacbioMapping

"""
    simple seeds selection for pathextends.
"""


def findseeds(graph):
    seeds = []
    for _, value in graph.vs.iteritems():
        seeds.append(value)
    seeds.sort(key=lambda x: x.length, reverse=True)
    return seeds


print "before load"
g = loadgraph(experimental.illumina_graph_file,
              experimental.illumina_contig_file,
              experimental.contig_refmap_file)
mseeds = findseeds(g)
start = mseeds[0]
print start

print "after load"
vertices = []
for _, value in g.vs.iteritems():
    vertices.append(value)
print "number of vertices", len(vertices)
# sorted(vertices, key=lambda v: min(v.positions_on_ref))
m = PacbioMapping(experimental.contig_pacbiomap_file)
print experimental.contig_pacbiomap_file
print len(m.contigToRead[str(start.vid)])


#now we sort the contigs length

def extendforward(graph, seed):
    current = seed
    print "extending"
    while True:
        neighbors = current.outv()
        if len(neighbors) > 0:
            current = neighbors[0]
            print current
        else:
            break

"""
ranking neighbours when genome is known
For development purpose only
"""


def rankneighbours_genome(g, current_vertex, neighbors):
    current_vertex.positions_on_ref
    overlap_vertices = {}
#    print "from:", current_vertex.positions_on_ref
    if len(current_vertex.positions_on_ref) == 0:
        print "no map", current_vertex
    for e in neighbors:
        overlap_vertices[e.v2] = e.ovl
    #    print "neighbors", [current_vertex.positions_on_ref[i] - e.ovl
    #                        for i in
    #                        range(0, len(current_vertex.positions_on_ref))]


def extendbackward(graph, seed):
    current = seed
    print "extending"
    while True:
        neighbors = current.innv()
        if len(neighbors) == 1:
            current = neighbors[0]
        elif len(neighbors) > 1:
            rankneighbours_genome(graph, current)
            current = neighbors[random.randint(0, len(neighbors) - 1)]
            print "move", current
            print len(neighbors)
        else:
            break

print "seed"
seeds = findseeds(g)
print "seed"
number_nomap = 0
for seed in seeds:
    if len(seed.positions_on_ref) == 0:
        number_nomap += 1
print "number of no mapping ", number_nomap
seed = seeds[0]
neighbors = seed.oute
if len(neighbors) > 0:
    rankneighbours_genome(g, seed, seed.oute)
#find peak number
