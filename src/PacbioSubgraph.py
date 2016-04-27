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

from PacbioAlignment import PacbioAlignment
from illumina_graph import Graph
from illumina_path import IlluminaPath
from utils import toggleStrand
import numpy as np

class PacbioSubgraph(Graph):

    def __init__(self, pacbio_id, illumina_graph,
                 pacbio_mapping, map_margin=0.5,
                 augment=False, offset_margin=0.1, min_vertex_length=0):
        Graph.__init__(self)
        self.pacbio_id = pacbio_id
        self.vertex_mapping = {}  # pacbio mapping corresponding to vertex
        self.map_offset = {}  # relative offset of each edge w.r.t. pacbio
        self.augmented_path = {}
        self.readArray = self.read_stack(pacbio_id, pacbio_mapping)
        ig = illumina_graph
        qStart = 0
        qEnd = 0
        qLength = 0
        for n in range(len(pacbio_mapping.readToContig[pacbio_id])):
            qStart = pacbio_mapping.readToContig[pacbio_id][n].qStart
            qEnd = pacbio_mapping.readToContig[pacbio_id][n].qEnd
            self.mark_coords(n, qStart, qEnd)
        self.covArray = self.readArray.sum(axis=0)
        
        #OLD CODE
        for i in range(len(pacbio_mapping.readToContig[pacbio_id])):
            for j in range(len(pacbio_mapping.readToContig[pacbio_id])):
                if j == i:
                    continue
                a1 = PacbioAlignment(pacbio_mapping.
                                     readToContig[pacbio_id][i])
                a2 = PacbioAlignment(pacbio_mapping.
                                     readToContig[pacbio_id][j])
                if a1.seqLen < min_vertex_length:
                    continue
                if a2.seqLen < min_vertex_length:
                    continue
                if a2.queryStartPos < a1.queryStartPos:
                    continue
                if a1.queryStrand == '+':
                    s1 = a1.seqStrand
                else:
                    s1 = toggleStrand(a1.seqStrand)
                if a2.queryStrand == '+':
                    s2 = a2.seqStrand
                else:
                    s2 = toggleStrand(a2.seqStrand)
                if ig.get_vid(a1.seqID, s1) not in self.vs.keys():
                    continue
                if ig.get_vid(a2.seqID, s2) not in self.vs.keys():
                    continue

                for e in ig.vs[ig.get_vid(a1.seqID, s1)].oute:
                    if e.v2.vid == ig.get_vid(a2.seqID, s2):
                        if a1.spansEdge(a2, e, offset_margin):
                            ne = self.add_edge(e.v1.vid, e.v2.vid, e.ovl)
                            ec = self.add_edge(e.v2.conj.vid,
                                               e.v1.conj.vid, e.ovl)
                            self.map_offset[ne] = a1.relative_offset(a2, ne.v1,
                                                                     ne.v2)
                            self.map_offset[ec] = a2.relative_offset(a1, ec.v1,
                                                                     ec.v2)
                            dbg_numedges += 1
                            break
        #print "[" + pacbio_id + "] Num edges: ", dbg_numedges
        if augment:
            self.augment(illumina_graph, offset_margin)
    

    #NEW CODE (04/27/16 DB)
    def read_stack(self, pacbio_id, pacbio_mapping):
        stack_size = len(pacbio_mapping.readToContig[pacbio_id])
        read_size = pacbio_mapping.readToContig[pacbio_id][0].qLength
        a = np.zeros((stack_size, read_size))
        return a
    def mark_coords(self, n, qStart, qEnd):
        for i in range(qStart, qEnd):
            self.readArray[n][i] = 1





    #OLD CODE
    def augment(self, illumina_graph, offset_margin):
        dbg_numpaths = 0
        for i in range(len(self.vs)):
            for j in range(len(self.vs)):
                if i == j:
                    continue
                vid1 = self.vs.keys()[i]
                vid2 = self.vs.keys()[j]
                if self.vs[vid2] in self.vs[vid1].outv():
                    continue
                a1 = self.vertex_mapping[vid1]
                a2 = self.vertex_mapping[vid2]
                offset = a1.relative_offset(a2,
                                            self.vs[vid1], self.vs[vid2])
                offsetc = a2.relative_offset(a1, self.vs[vid2].conj,
                                             self.vs[vid1].conj)
                if offset is False or offset < 0:
                    continue
                if offsetc is False or offsetc < 0:
                    continue
                min_offset = (1 - offset_margin) * offset
                max_offset = (1 + offset_margin) * offset
                if self.find_path(vid1, vid2, min_offset, max_offset):
                    continue
                p = illumina_graph.find_path(vid1, vid2,
                                             min_offset, max_offset)
                if p is not None and len(p) > 1:
                    dbg_numpaths += 1
                    ip = IlluminaPath(p)
                    if a1.spansEdge(a2, ip, offset_margin):
                        continue
                    ne = self.add_edge(ip.v1.vid, ip.v2.vid, ip.ovl)
                    ec = self.add_edge(ip.v2.conj.vid, ip.v1.conj.vid,
                                       ip.ovl)
                    self.map_offset[ne] = a1.relative_offset(a2, ne.v1, ne.v2)
                    self.map_offset[ec] = a2.relative_offset(a1, ec.v1, ec.v2)
                    self.augmented_path[ne] = ip
                    self.augmented_path[ec] = ip.conj()
        #print "Num paths: ", dbg_numpaths
        if dbg_numpaths > 0:
            return True
        return False
