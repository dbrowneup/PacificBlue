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

from illumina_graph import Graph
from PacbioSubgraph import PacbioSubgraph
from PacbioAlignment import PacbioAlignment
from utils import toggleStrand, parmap

from sets import Set
import math


class LongContigGraph(Graph):

    def __init__(self, contigGraph, mapping, map_margin=0.7,
                 augment_paths=True, offset_margin=0.1, num_threads=1,
                 read_thresholds=(3,10)):
        Graph.__init__(self)
        self.mapping = mapping
        self.contig_graph = contigGraph
        self.read_count = {}
        self.map_margin = map_margin
        self.offset_margin = offset_margin
        self.num_reads = 0
        self.num_edged_graphs = 0
        self.num_multivertices = 0
        self.num_singlevertices = 0
        self.num_threads = num_threads
        self.read_count_threshold_high = read_thresholds[1]
        self.read_count_threshold_low = read_thresholds[0]
        vertex_length_schedule = [100000, 50000, 25000, 12500, 6250, 1000]
        previous_vertex_length = 10000000000

        mapping_reads = self.mapping.readToContig.keys()
        batch_size = int(math.ceil(float(len(mapping_reads)) /
                               (self.num_threads))) + 1
        print "Number of threads = ", self.num_threads
        print "Batch size: ", batch_size
        for min_vertex_length in vertex_length_schedule:
            print "Contig length threshold = ", min_vertex_length

            def parallel_subgraphs(thread_id=0):
                num_reads = 0
                num_edged_graphs = 0
                num_multivertices = 0
                num_singlevertices = 0
                vertex_list = []
                edge_list = []
                for r in mapping_reads[thread_id * batch_size:
                                       min((thread_id + 1) * batch_size,
                                           len(mapping_reads))]:
                    if len(mapping.readToContig[r]) < 2:
                        continue
                    sg = PacbioSubgraph(r, contigGraph, self.mapping,
                                        map_margin, True, offset_margin,
                                        min_vertex_length)
                    num_reads += 1
                    if len(sg.es) > 0:
                        num_edged_graphs += 1
                    if len(sg.vs) > 2:
                        num_multivertices += 1
                    if len(sg.vs) > 0:
                        num_singlevertices += 1
                    for v1id in sg.vs:
                        if sg.vertex_mapping[v1id].seqLen < min_vertex_length:
                            continue
                        v1 = sg.vs[v1id]
                        for e in sg.vs[v1id].oute:
                            v2 = e.v2
                            v2id = v2.vid
                            if (sg.vertex_mapping[v2id].seqLen <
                                min_vertex_length):
                                continue
                            if not sg.vertex_mapping[v1id].\
                            ConcordantDirection(sg.vertex_mapping[v2id],
                                                v1, v2):
                                continue
                            reloff = sg.vertex_mapping[v1id].\
                            relative_offset(sg.vertex_mapping[v2id],
                                            sg.vs[v1id], sg.vs[v2id])
                            v1idc = v1.conj.vid
                            v2idc = v2.conj.vid
                            reloffc = sg.vertex_mapping[v2idc].\
                            relative_offset(sg.vertex_mapping[v1idc],
                                            sg.vs[v2idc], sg.vs[v1idc])
                            if reloff > 0 and reloffc > 0:
                                if v1id not in self.vs or v2id not in self.vs:
                                    vertex_list.append((v1id, v1.conj.vid,
                                                        v1.length, v1.cov))
                                    vertex_list.append((v1.conj.vid, v1id,
                                                        v1.length, v1.cov))
                                    vertex_list.append((v2id, v2.conj.vid,
                                                        v2.length, v2.cov))
                                    vertex_list.append((v2.conj.vid, v2id,
                                                        v2.length, v2.cov))
                                if (v1.vid in self.vs and
                                    v1.length >= previous_vertex_length):
                                    long_neighbors = False
                                    for ev1 in self.vs[v1.vid].oute:
                                        if (ev1.v2.length >
                                            previous_vertex_length):
                                            long_neighbors = True
                                            break
                                    if long_neighbors:
                                        continue
                                if (v2.vid in self.vs and
                                    v2.length >= previous_vertex_length):
                                    long_neighbors = False
                                    for ev2 in self.vs[v2.vid].inne:
                                        if (ev2.v1.length >
                                            previous_vertex_length):
                                            long_neighbors = True
                                            break
                                    if long_neighbors:
                                        continue
                                edge_list.append(((v1id, v2id, e.ovl),
                                                  (v2.conj.vid, v1.conj.vid,
                                                   e.ovl)))
                return (vertex_list, edge_list,
                        (num_reads, num_edged_graphs,
                         num_multivertices, num_singlevertices))

            print "Adding new vertices and edges to skeleton"
            new_graphs = parmap(parallel_subgraphs, range(self.num_threads))
            for vs, es, rstats in new_graphs:
                for vv in vs:
                    self.add_vertex(vv[0], vv[1], vv[2], vv[3])
                for ee in es:
                    ne = self.add_edge(ee[0][0], ee[0][1], ee[0][2])
                    ec = self.add_edge(ee[1][0], ee[1][1], ee[1][2])
                    if ne.eid in self.read_count:
                        self.read_count[ne.eid] += 1
                    else:
                        self.read_count[ne.eid] = 1
                    if ec.eid in self.read_count:
                        self.read_count[ec.eid] += 1
                    else:
                        self.read_count[ec.eid] = 1
                self.num_reads += rstats[0]
                self.num_edged_graphs += rstats[1]
                self.num_multivertices += rstats[2]
                self.num_singlevertices += rstats[3]

            edge_keys = self.es.keys()
            for eid in edge_keys:
                if eid not in self.es:
                    continue
                if self.read_count[eid] < self.read_count_threshold_low:
                    self.remove_edge(self.es[eid])
                    continue
                insignificant = -1
                not_best = -1
                if self.read_count[eid] >= self.read_count_threshold_high:
                    continue
                for e1 in self.es[eid].v1.oute:
                    if e1.eid == eid:
                        continue
                    if self.read_count[e1.eid] >= 2 * self.read_count[eid]:
                        insignificant += 1
                    if self.read_count[e1.eid] > self.read_count[eid]:
                        not_best += 1

                for e2 in self.es[eid].v2.inne:
                    if e2.eid == eid:
                        continue
                    if self.read_count[e2.eid] >= 2 * self.read_count[eid]:
                        insignificant += 1
                    if self.read_count[e2.eid] > self.read_count[eid]:
                        not_best += 1
                if insignificant >= 0 and not_best > 0:
                    self.remove_edge(self.es[eid])

            self.remove_transitive_edges(offset_margin)

            self.resolve_repeats()

            self.extend_ends(10, min_vertex_length)

            self.remove_transitive_edges(offset_margin)
            self.resolve_repeats()

            previous_vertex_length = min_vertex_length

    def extend_ends(self, read_support=5, min_vertex_length=-1):
        if min_vertex_length == -1:
            min_vertex_length = self.min_vertex_length
        end_vertices = Set([])
        for vid in self.vs:
            if len(self.vs[vid].oute) == 0 and len(self.vs[vid].inne) <= 1:
                end_vertices.add(self.vs[vid])
        for vid in self.contig_graph.vs:
            if (vid not in self.vs and
                self.contig_graph.vs[vid].length >= min_vertex_length):
                end_vertices.add(self.contig_graph.vs[vid])
        end_vertex_ids = Set([v.vid for v in end_vertices])

        branch_vertices = Set([])
        for vid in self.vs:
            if vid in end_vertex_ids:
                continue
            if len(self.vs[vid].oute) > 1 or len(self.vs[vid].inne) > 1:
                branch_vertices.add(self.vs[vid])
        branch_vertex_ids = Set([v.vid for v in branch_vertices])

        end_vertex_pairs = {}
        mixed_vertex_pairs = {}
        neighbors = {}

        for v in end_vertices:
            if v.Name() not in self.mapping.contigToRead:
                continue
            for m1 in self.mapping.contigToRead[v.Name()]:
                a1 = PacbioAlignment(m1, self.mapping.fileFormat)
                if a1.seqStrand == a1.queryStrand:
                    s1Strand = '+'
                else:
                    s1Strand = '-'
                if v.vid < 0:
                    if s1Strand == '+':
                        if a1.seqStartPos > a1.queryStartPos:
                            continue
                    else:
                        if a1.seqStartPos > a1.queryLen - a1.queryEndPos:
                            continue
                    s1Strand = toggleStrand(s1Strand)
                else:
                    if s1Strand == '+':
                        if a1.seqLen - a1.seqEndPos >\
                        a1.queryLen - a1.queryEndPos:
                            continue
                    else:
                        if a1.seqLen - a1.seqEndPos > a1.queryStartPos:
                            continue

                for m2 in self.mapping.readToContig[a1.queryID]:
                    a2 = PacbioAlignment(m2, self.mapping.fileFormat)
                    if a2.seqID == a1.seqID:
                        continue
                    if a2.seqStrand == a2.queryStrand:
                        s2Strand = '+'
                    else:
                        s2Strand = '-'
                    a1s = a1.queryStartPos
                    a1e = a1.queryEndPos
                    a2s = a2.queryStartPos
                    a2e = a2.queryEndPos
                    if ((a1s - a2s) * (a1e - a2s) < 0 and
                        (a1s - a2e) * (a1e - a2e)) < 0:
                        continue
                    if ((a2s - a1s) * (a2e - a1s) < 0 and
                        (a2s - a1e) * (a2e - a1e)) < 0:
                        continue
                    if s1Strand != s2Strand:
                        if s2Strand == '+':
                            if a2.seqLen - a2.seqEndPos >\
                            a2.queryLen - a2.queryEndPos:
                                continue
                        else:
                            if a2.seqLen - a2.seqEndPos > a2.queryStartPos:
                                continue
                        reloff1 = a1.relative_offset(a2, v, self.contig_graph.
                                                     vs[int(a2.seqID)].conj)
                        reloff2 = a2.relative_offset(a1, self.contig_graph.
                                                     vs[int(a2.seqID)], v.conj)
                        if reloff1 < 0 or reloff2 < 0:
                            continue
                        if reloff1 < v.length:
                            if v.vid > 0:
                                if reloff1 < a1.seqStartPos:
                                    continue
                            else:
                                if v.length - a1.seqEndPos > reloff1:
                                    continue
                            if a2.seqStartPos > reloff2:
                                continue

                        if int(a2.seqID) in end_vertex_ids:
                            if (v.vid, int(a2.seqID)) in end_vertex_pairs:
                                end_vertex_pairs[(v.vid, int(a2.seqID))].\
                                append((a1, a2))
                            else:
                                end_vertex_pairs[(v.vid, int(a2.seqID))] =\
                                [(a1, a2)]
                        elif int(a2.seqID) in self.vs:
                            if (v.vid, int(a2.seqID)) in mixed_vertex_pairs:
                                mixed_vertex_pairs[(v.vid, int(a2.seqID))].\
                                append((a1, a2))
                            else:
                                mixed_vertex_pairs[(v.vid, int(a2.seqID))] =\
                                [(a1, a2)]
                    if s1Strand == s2Strand:
                        if s2Strand == '+':
                            if a2.seqStartPos > a2.queryStartPos:
                                continue
                        else:
                            if (a2.seqStartPos >
                                a2.queryLen - a2.queryEndPos):
                                continue
                        v2id = -1 * int(a2.seqID)
                        reloff1 = a1.relative_offset(a2, v, self.contig_graph.
                                                     vs[v2id].conj)
                        reloff2 = a2.relative_offset(a1, self.contig_graph.
                                                     vs[v2id], v.conj)
                        if reloff1 < 0 or reloff2 < 0:
                            continue
                        if reloff1 < v.length:
                            if v.vid > 0:
                                if reloff1 < a1.seqStartPos:
                                    continue
                            else:
                                if v.length - a1.seqEndPos > reloff1:
                                    continue
                            if a2.seqLen - a2.seqEndPos > reloff2:
                                continue
                        if v2id in end_vertex_ids:
                            if (v.vid, v2id) in end_vertex_pairs:
                                end_vertex_pairs[(v.vid, v2id)].\
                                append((a1, a2))
                            else:
                                end_vertex_pairs[(v.vid, v2id)] = [(a1, a2)]
                        elif v2id in self.vs:
                            if (v.vid, v2id) in mixed_vertex_pairs:
                                mixed_vertex_pairs[(v.vid, v2id)].\
                                append((a1, a2))
                            else:
                                mixed_vertex_pairs[(v.vid, v2id)] =\
                                [(a1, a2)]

        print '# putative vertex pairs for gap filling:', len(end_vertex_pairs)

        for vp in end_vertex_pairs.keys() + mixed_vertex_pairs.keys():
            if vp[0] in neighbors:
                neighbors[vp[0]].add(vp[1])
            else:
                neighbors[vp[0]] = Set([vp[1]])

        best_match = {}
        new_pairs = {}
        end_matches = {}
        for v1id in neighbors:
            end_neighbors = [v2id for v2id in neighbors[v1id]
                             if v2id in end_vertex_ids]
            end_matches[v1id] = [(v2id, len(end_vertex_pairs[(v1id, v2id)]))
                                 for v2id in neighbors[v1id]
                                 if v2id in end_vertex_ids]
            end_matches[v1id].sort(lambda x, y: x[1] - y[1])

            max1 = [0, 0]
            max2 = [0, 0]
            for v2id in end_neighbors:
                if len(end_vertex_pairs[(v1id, v2id)]) > max1[1]:
                    max2 = max1
                    max1 = (v2id, len(end_vertex_pairs[(v1id, v2id)]))
                elif len(end_vertex_pairs[(v1id, v2id)]) > max2[1]:
                    max2 = (v2id, len(end_vertex_pairs[(v1id, v2id)]))
            if (max1[1] > self.read_count_threshold_low and
                max2[1] < self.read_count_threshold_high and
                max1[1] >= max2[1] + max(max2[1],
                                         self.read_count_threshold_low)):
                best_match[v1id] = max1
            if max1[1] > self.read_count_threshold_high:
                continue
            branch_neighbors = [v2id for v2id in neighbors[v1id]
                                if v2id in branch_vertex_ids]
            branchmax1 = [0, 0]
            branchmax2 = [0, 0]
            for v2id in branch_neighbors:
                if len(mixed_vertex_pairs[(v1id, v2id)]) > branchmax1[1]:
                    branchmax2 = branchmax1
                    branchmax1 = (v2id, len(mixed_vertex_pairs[(v1id, v2id)]))
                elif len(mixed_vertex_pairs[(v1id, v2id)]) > branchmax2[1]:
                    branchmax2 = (v2id, len(mixed_vertex_pairs[(v1id, v2id)]))
            if branchmax2[1] <= max1[1]:
                nextbest = max1
            else:
                nextbest = branchmax2
            if (branchmax1[1] > self.read_count_threshold_low and
                nextbest[1] < self.read_count_threshold_high and
                branchmax1[1] >= (nextbest[1] +
                                  max(nextbest[1],
                                      self.read_count_threshold_low))):
                best_match[v1id] = branchmax1

        for v1id in best_match:
            v2id = best_match[v1id][0]
            if v2id in end_vertex_ids:
                if best_match[v1id][1] == end_matches[v2id][-1][1]:
                    minvid = min(v1id, v2id)
                    maxvid = max(v1id, v2id)
                    new_pairs[(minvid, maxvid)] = best_match[v1id][1]
                    continue
            if v2id in branch_vertex_ids:
                new_pairs[(v1id, v2id)] = best_match[v1id][1]
                continue

        num_gaps_filled = 0
        for vp in new_pairs:
            v0 = self.contig_graph.vs[vp[0]]
            v1 = self.contig_graph.vs[vp[1]]
            if vp[0] not in self.vs:
                self.add_vertex(v0.vid, v0.conj.vid, v0.length, v0.cov)
                self.add_vertex(v0.conj.vid, v0.vid, v0.length, v0.cov)
            if vp[1] not in self.vs:
                self.add_vertex(v1.vid, v1.conj.vid, v1.length, v1.cov)
                self.add_vertex(v1.conj.vid, v1.vid, v1.length, v1.cov)
            if vp in end_vertex_pairs:
                vpaligns = end_vertex_pairs[vp]
            else:
                vpaligns = mixed_vertex_pairs[vp]
            reloffs = [ap[0].relative_offset(ap[1], v0, v1.conj)
                       for ap in vpaligns]
            ne = self.add_edge(v0.vid, v1.conj.vid,
                          reloffs[len(reloffs) / 2] - v0.length)
            ec = self.add_edge(v1.vid, v0.conj.vid,
                          reloffs[len(reloffs) / 2] - v0.length)
            num_gaps_filled += 1
            self.read_count[ec.eid] = len(vpaligns)
            self.read_count[ne.eid] = len(vpaligns)
        print "# gaps filled = ", num_gaps_filled

    def resolve_repeats(self):
        vskeys = [vid for vid in self.vs.keys() if vid > 0]
        vskeys.sort(lambda x, y: self.vs[x].length - self.vs[y].length)
        num_spanned = 0
        for vid in vskeys:
            v = self.vs[vid]
            if len(v.inne) <= 1 and len(v.oute) <= 1:
                continue
            if len(v.inne) < 1 or len(v.oute) < 1:
                continue
            supported_pairs = {}
            for e1 in v.inne:
                for e2 in v.oute:
                    v1 = e1.v1
                    v2 = e2.v2
                    offset = v1.length + e1.ovl + v.length + e2.ovl
                    ovl = e1.ovl + v.length + e2.ovl
                    min_offset = (1 - self.offset_margin) * offset
                    max_offset = (1 + self.offset_margin) * offset

                    for m1 in self.mapping.contigToRead[v1.Name()]:
                        a1 = PacbioAlignment(m1, self.mapping.fileFormat)
                        if not a1.longAlignment(self.map_margin):
                            continue
                        if len(self.mapping.readToContig[a1.queryID]) < 2:
                            continue
                        a2 = None
                        for m2 in self.mapping.readToContig[a1.queryID]:
                            a2 = PacbioAlignment(m2, self.mapping.fileFormat)
                            if a2.seqID != v2.Name():
                                continue
                            if a2.longAlignment(self.map_margin):
                                break
                        if a2 is None or a2.seqID != v2.Name():
                            continue
                        reloff = a1.relative_offset(a2, v1, v2)
                        if reloff == False:
                            continue
                        if reloff >= min_offset and reloff <= max_offset:
                            spanning1 = False
                            spanning2 = False
                            if v1.Strand() == '+':
                                if a1.seqStartPos < v1.length + e1.ovl:
                                    spanning1 = True
                            if v1.Strand() == '-':
                                if a1.seqEndPos > -1 * e1.ovl:
                                    spanning1 = True
                            if v2.Strand() == '+':
                                if a2.seqEndPos > -1 * e2.ovl:
                                    spanning2 = True
                            if v2.Strand() == '-':
                                if a2.seqStartPos < v2.length + e2.ovl:
                                    spanning2 = True
                            if spanning1 and spanning2:
                                if (v1.vid, v2.vid) in supported_pairs:
                                    supported_pairs[(v1.vid, v2.vid)][0] += 1
                                else:
                                    supported_pairs[(v1.vid, v2.vid)] =\
                                    [1, ovl]
            for vp in supported_pairs:
                if (len(self.vs[vp[0]].inne) > 1 or
                    len(self.vs[vp[0]].oute) > 1):
                    continue
                if (len(self.vs[vp[1]].inne) > 1 or
                    len(self.vs[vp[1]].oute) > 1):
                    continue
                if supported_pairs[vp][0] >= self.read_count_threshold_high:
                    self.add_transitive_edge(vp, v, supported_pairs[vp])
                    continue
                if supported_pairs[vp][0] < self.read_count_threshold_low:
                    continue
                v1max2 = 0
                v2max2 = 0
                for vp2 in supported_pairs:
                    if vp2 == vp:
                        continue
                    if vp2[0] == vp[0]:
                        if supported_pairs[vp2][0] > v1max2:
                            v1max2 = supported_pairs[vp2][0]
                    if vp2[1] == vp[1]:
                        if supported_pairs[vp2][0] > v2max2:
                            v2max2 = supported_pairs[vp2][0]
                if ((v1max2 + max(v1max2, self.read_count_threshold_low) <=
                     supported_pairs[vp][0]) and
                    (v2max2 + max(v1max2, self.read_count_threshold_low) <=
                     supported_pairs[vp][0])):
                    self.add_transitive_edge(vp, v, supported_pairs[vp])
                    num_spanned += 1
        print "# branches resolved = ", num_spanned

    def add_transitive_edge(self, vp, v, pair_information):
        read_count = pair_information[0]
        ovl = pair_information[1]
        for e1 in self.vs[vp[0]].oute:
            if e1.v2.vid == v.vid:
                self.remove_edge(e1)
        for e2 in self.vs[vp[1]].inne:
            if e2.v1.vid == v.vid:
                self.remove_edge(e2)
        ne = self.add_edge(vp[0], vp[1], ovl)
        ec = self.add_edge(self.vs[vp[1]].conj.vid,
                           self.vs[vp[0]].conj.vid, ovl)
        self.read_count[ne.eid] = read_count
        self.read_count[ec.eid] = read_count
        return ne
