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

# from illumina_graph import Graph, OEdge
# import sys


class IlluminaPath():

    def __init__(self, path):
#         self.graph = Graph()
        self.path = path
        self.v1 = path[0].v1
        self.v2 = path[-1].v2
        self.ovl = -1 * path[0].v1.length
        for e in path:
            self.ovl += e.v1.length + e.ovl

#     def add_vertex(self, v, v_conj):
#         if len(self.graph.vs) > 0:
#             return False
#         assert (v.conj == v_conj) and (v == v_conj.conj),\
#             "non-cojugate vertices"
#         self.graph.vs[v.vid] = v
#         self.graph.vs[v_conj.vid] = v_conj
#         self.source = v
#         self.sink = v
#         self.ovl = -1 * v.length

#     def extend(self, e):
#         if len(self.graph.vs) == 0:
#             v1 = self.graph.add_vertex(e.v1.vid, e.v1.conj.vid,
#                                   e.v1.length, e.v1.cov)
#             v1conj = self.graph.add_vertex(e.v1.conj.vid, e.v1.vid,
#                                   e.v1.length, e.v1.cov)
#             self.ovl = -1 * e.v1.length
#             self.source = v1
#             self.sink = v1
#         if self.graph.has_edge(e.v1.vid, e.v2.vid, e.ovl):
#             return True
#         if (e.v1.vid == self.sink.vid):
#             v2 = self.graph.add_vertex(e.v2.vid, e.v2.conj.vid,
#                                   e.v2.length, e.v2.cov)
#             self.graph.add_vertex(e.v2.conj.vid, e.v2.vid,
#                                   e.v2.length, e.v2.cov)
#             v1 = self.graph.vs[e.v1.vid]
#             v2 = self.graph.vs[e.v2.vid]
#             self.graph.add_edge(v1.vid, v2.vid, e.ovl, e.eid)
#             if e.conj is not None:
#                 self.graph.add_edge(v2.conj.vid, v1.conj.vid,
#                                     e.conj.ovl, e.conj.eid)
#             else:
#                 sys.exit("Error: conjugate edge absent")
#             self.sink = v2
#             self.ovl = self.ovl + e.v1.length + e.ovl
#             return True
#         if (e.v2 == self.source):
#             v1 = self.graph.add_vertex(e.v1.vid, e.v1.conj.vid,
#                                   e.v1.length, e.v1.cov)
#             self.graph.add_vertex(e.v1.conj.vid, e.v1.vid,
#                                   e.v1.length, e.v1.cov)
#             v1 = self.graph.vs[e.v1.vid]
#             v2 = self.graph.vs[e.v2.vid]
#             self.graph.add_edge(v1.vid, v2.vid, e.ovl, e.eid)
#             if e.conj is not None:
#                 self.graph.add_edge(v2.conj.vid, v1.conj.vid,
#                                     e.conj.ovl, e.conj.eid)
#             else:
#                 sys.exit("Error: conjugate edge absent")
#             self.source = v1
#             self.ovl = self.ovl + e.v2.length + e.ovl
#             return True
#         if (e.v1.conj == self.source):
#             v2 = self.graph.add_vertex(e.v2.vid, e.v2.conj.vid,
#                                   e.v2.length, e.v2.cov)
#             self.graph.add_vertex(e.v2.conj.vid, e.v2.vid,
#                                   e.v2.length, e.v2.cov)
#             v1 = self.graph.vs[e.v1.vid]
#             v2 = self.graph.vs[e.v2.vid]
#             self.graph.add_edge(v1.vid, v2.vid, e.ovl, e.eid)
#             if e.conj is not None:
#                 self.graph.add_edge(v2.conj.vid, v1.conj.vid,
#                                     e.conj.ovl, e.conj.eid)
#             else:
#                 sys.exit("Error: conjugate edge absent")
#             self.source = v2.conj
#             self.ovl = self.ovl + e.v1.length + e.ovl
#             return True
#         if (e.v2.conj == self.sink):
#             v1 = self.graph.add_vertex(e.v1.vid, e.v1.conj.vid,
#                                   e.v1.length, e.v1.cov)
#             self.graph.add_vertex(e.v1.conj.vid, e.v1.vid,
#                                   e.v1.length, e.v1.cov)
#             v1 = self.graph.vs[e.v1.vid]
#             v2 = self.graph.vs[e.v2.vid]
#             self.graph.add_edge(v1.vid, v2.vid, e.ovl, e.eid)
#             if e.conj is not None:
#                 self.graph.add_edge(v2.conj.vid, v1.conj.vid,
#                                     e.conj.ovl, e.conj.eid)
#             else:
#                 sys.exit("Error: conjugate edge absent")
#             self.sink = v1.conj
#             self.ovl = self.ovl + e.v2.length + e.ovl
#             return True
#         return False

    def conj(self):

        return IlluminaPath([e.conj for e in self.path[::-1]])

#         ipconj = IlluminaPath()
#         ipconj.graph = self.graph
#         ipconj.source = self.sink.conj
#         ipconj.sink = self.source.conj
#         ipconj.ovl = self.ovl
#         return ipconj
