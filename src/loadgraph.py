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

#pathextend for pacbio contigs
from illumina_graph import Graph
from PacbioAlignment import PacbioAlignment
from PacbioMapping import PacbioMapping
#class PathExtend:
#    def __init(self, pacbiomap, graph):
#        self.readmap = pacbiomap
#        self.graph = graph


def loadgraph(illumina_graph_file, illumina_contig_file, contig_refmap_file):

    vertices = []

    m_ref = PacbioMapping(contig_refmap_file, 'blastn')
    g = Graph()
    g.load(illumina_graph_file, illumina_contig_file)

    for key, _ in g.vs.iteritems():
        vertices.append(key)

    for key, value in m_ref.readToContig.iteritems():
        positions = []
        positions_rc = []
        if len(value) > 0:
            for i in range(len(value)):
                pa = PacbioAlignment(value[i].strip(), 'blastn')
                position = pa.seqStartPos
                if pa.seqStrand == '-':
                    position = -1 * position

                positions.append(position)
            if abs(int(key)) == 657:
                print 'position', position

            g.add_positions_ref(int(key), positions)
#            if int(key) == 657:
#                print "hello world", positions
#            print "add", key
            length_contig = g.vs[int(key)].length
            positions_rc = [(abs(pos) - length_contig) * -1 * cmp(pos, 0)
                            for pos in positions]
            g.add_positions_ref(-1 * int(key), positions_rc)

#            print "length:", length_contig
#            print "adding position", key

#    print "list", len(vertices)
#    for id in vertices:
#        if id not in inmapping and -1*id not in inmapping:
#            print "Length", id, g.vs[id].length
    return g
