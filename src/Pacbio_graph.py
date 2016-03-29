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

from abstract_graph import Abstract_Graph, Abstract_Edge, Abstract_Vertex
#import utils
import logging
import parser
import sys
from utils import conjugate


class OVertex(Abstract_Vertex):
    def __init__(self, vid, conj, vlen, coverage):
        Abstract_Vertex.__init__(self, vid)
        conjugate(self, conj)
        self.length = vlen
        self.seq = None
        self.cov = coverage
        self.positions_on_ref = []

    def __repr__(self):
        return "V: %d(l: %d)(c: %d)" % (self.vid, self.length, self.cov)

    def Strand(self):
        if self.vid > 0:
            return '+'
        if self.vid < 0:
            return '-'

    def Name(self):
        return str(abs(self.vid))


class OEdge(Abstract_Edge):
    def __init__(self, edgeid, v1, v2, overlap):
        Abstract_Edge.__init__(self, v1, v2, edgeid)
        self.ovl = overlap
        self.conj = None

    def __repr__(self):
        return "E%d(%d)--> (%d),%d " % \
            (self.eid, self.v1.vid, self.v2.vid, self.ovl)


class Graph(Abstract_Graph):
    def __init__(self):
        Abstract_Graph.__init__(self)
        self.max_eid = 0
        self.logger = logging.getLogger('Algae')

    def add_vertex(self, vid, conj_id, vlen, cov):
        conj = self.vs.get(conj_id, None)
        v = OVertex(vid, conj, vlen, cov)
        self.vs[vid] = v
        return v

#     def include_vertex(self, v, v_conj):
#         if v.vid in self.vs:
#             assert self.vs[v.vid] == v, "Vertex duplicate"
#         if v_conj.vid in self.vs:
#             assert self.vs[v_conj.vid] == v_conj, "Vertex duplicate"
#         self.vs[v.vid] = v
#         self.vs[v_conj.vid] = v_conj
#         return v

    def add_edge(self, v1id, v2id, overlap, eid=None):
        v1 = self.vs[v1id]
        v2 = self.vs[v2id]
        if v2 in v1.outv():
            for e in v1.oute:
                if e.v2 == v2 and (eid is None or e.eid == eid):
                    return e
            sys.exit('Edge not found!')

        if eid is None:
            while self.max_eid in self.es.keys():
                self.max_eid += 1
            neweid = self.max_eid
            self.max_eid += 1
        else:
            if eid in self.es.keys():
                sys.exit("EID exists")
            neweid = eid
        e = OEdge(neweid, v1, v2, overlap)
        self.es[e.eid] = e
        #if the conjugate already exists, add conjugate
        v1conj = v1.conj
        v2conj = v2.conj
        if v1conj in v2conj.outv():
            ej = None
            for ec in v2conj.oute:
                if ec.v2 == v1conj:
                    ej = ec
            if not ej:
                sys.exit('Conjugate edge not found!')
            conjugate(e, ej)
        return e

    def has_edge(self, v1id, v2id, overlap, offset_margin=0):
        if v1id not in self.vs or v2id not in self.vs:
            return False
        if v2id not in self.vs[v1id].outv():
            return False
        for e in self.vs[v1id].oute:
            if e.v2.vid == v2id and (overlap >= (1 - offset_margin) * e.ovl and
                                     overlap <= (1 + offset_margin) * e.ovl):
                return True
        return False

    def find_path(self, v1id, v2id, min_offset, max_offset, depth=0):
        if max_offset < 0:
            return None
        if depth > 10:
            return None
        v1 = self.vs[v1id]
        if v1id == v2id:
            if min_offset <= 0 and max_offset >= 0:
                return []
        for ve in v1.oute:
            off = v1.length + ve.ovl
            p = self.find_path(ve.v2.vid, v2id, min_offset - off,
                               max_offset - off, depth + 1)
            if p is not None:
                return [ve] + p
        return None

    def add_seq(self, vid, seq):
        self.vs[vid].seq = seq

    def load_dna(self, filename):
        for contigid, seq in parser.contigs_sequence(filename):
            self.add_seq(contigid, seq)
#            print contigid , seq

    def add_cov(self, vid, cvr):
        self.vs[vid].cvr = cvr

    def path_seq(self, path):
        vertices = []
        overlapping = []
        for e in path:
            vertices.append(e.v1)
            overlapping.append(e.ovl)

        vertices.append(path[-1].v2)
        seq = vertices[0].seq
        #do not consider gaps
        for i in range(1, len(vertices)):
            seq += vertices[i].seq[-1 * overlapping[i - 1]:]
        return seq

    def load(self, graph_filename, seq_filename):
        for vid, conj, length, cov in parser.graph_vertices(graph_filename):
            self.add_vertex(vid, conj, length, cov)
            self.add_vertex(conj, vid, length, cov)
#            print va, vb, va.conj, vb.conj
        for v1, v2, overlap in parser.graph_edges(graph_filename):
            self.add_edge(v1, v2, overlap)
#            print overlap

    def get_vid(self, name, strand):
        if strand == '+':
            return int(name)
        else:
            return -1 * int(name)
    """
    vid: vertex id
    positions:  if positions < 0 it means in the reverse complement

    """
    def add_positions_ref(self, vid, positions):
        self.vs[vid].positions_on_ref = positions[:]
        print "ok", self.vs[vid].positions_on_ref

g = Graph()
#g.load("data/spalgae-contigs.dot",None)
#branches = 0
#small_branches = 0
#
#for i in g.vs.values():
#    if len(i.outv()) > 1 and len(i.innv())> 1:
#
##        print "multiple"
#        print i.length
#        if i.length < 100:
#            small_branches +=1
#        branches +=1
##        print "Middle", i
##        print "Out" , i.outv()
##        print "Inn", i.innv()
#portion = "small contigs", small_branches*1.0/branches
#print "small contigs portion", portion

g.load_dna("data/secoli_pacbio_contigs_mapping.fasta.m4")
print g.path_seq([g.es[0]])

for i in g.vs.values():
    print i

for i in g.es.values():
    print i
