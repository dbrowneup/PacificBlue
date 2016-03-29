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

from PacbioMapping import PacbioMapping
from PacbioAlignment import PacbioAlignment
from sets import Set


class MappingValidation:

    def __init__(self, c2rFileName, p2cFileName, p2rFileName,
                 c2rFileFormat='blastn', p2cFileFormat='m4',
                 p2rFileFormat='m4'):
        self.c2rMapping = PacbioMapping(c2rFileName, c2rFileFormat, 0.9)
        self.p2rMapping = PacbioMapping(p2rFileName, p2rFileFormat)
        self.p2cMapping = PacbioMapping(p2cFileName, p2cFileFormat)
        self.c2rFileFormat = c2rFileFormat
        self.p2rFileFormat = p2rFileFormat
        self.p2cFileFormat = p2cFileFormat
        self.alphabet = {}

        def last_pos_cmp(m1, m2):
            a1 = PacbioAlignment(m1, self.c2rFileFormat)
            a2 = PacbioAlignment(m2, self.c2rFileFormat)
            if a1.seqEndPos != a2.seqEndPos:
                return a1.seqEndPos - a2.seqEndPos
            if a1.seqStartPos != a2.seqStartPos:
                return a1.seqStartPos - a2.seqStartPos
            return a1.queryStartPos - a2.queryStartPos

        for c in self.c2rMapping.contigToRead:
            self.alphabet[c] = self.c2rMapping.contigToRead[c][:]
            self.alphabet[c].sort(last_pos_cmp)

    def valid_edge(self, v1, v2, offset, offset_margin=0.1):
        # TODO
        for m1 in self.mapping.readToContig[v1.Name()]:
            for m2 in self.mapping.readToContig[v2.Name()]:
                a1 = PacbioAlignment(m1, self.fileFormat)
                a2 = PacbioAlignment(m2, self.fileFormat)
                qoff = a1.query_offset(a2, v1, v2)
                if qoff < (1 + offset_margin) * offset:
                    if qoff > (1 - offset_margin) * offset:
                        return True
        return False

    def valid_trio(self, p2c, c2r, p2r, offset_margin=0.1):
        strand = 1
        if p2c.queryStrand != p2c.seqStrand:
            strand = strand * -1
        if c2r.queryStrand != c2r.seqStrand:
            strand = strand * -1
        if p2r.queryStrand == p2r.seqStrand and strand == -1:
            return False
        elif p2r.queryStrand != p2r.seqStrand and strand == 1:
            return False

        if strand == 1:
            if p2c.queryStrand == p2c.seqStrand:
                cpos = p2c.seqStartPos - p2c.queryStartPos
            else:
                cpos = p2c.seqEndPos + p2c.queryStartPos
        else:
            if p2c.queryStrand == p2c.seqStrand:
                cpos = p2c.seqEndPos + (p2c.queryLen - p2c.queryEndPos)
            else:
                cpos = p2c.seqStartPos - (p2c.queryLen - p2c.queryEndPos)

        if c2r.queryStrand == c2r.seqStrand:
            rpos = c2r.seqStartPos - c2r.queryStartPos + cpos
        else:
            rpos = c2r.seqEndPos - (cpos - c2r.queryStartPos)

        fpos = 0
        if p2r.seqStrand == p2r.queryStrand:
            fpos = p2r.seqStartPos - p2r.queryStartPos
        else:
            fpos = p2r.seqStartPos - (p2r.queryLen - p2r.queryEndPos)

        return (abs(fpos - rpos) < offset_margin * p2r.queryLen)

    def trio_stats(self, map_margin=0.5):
        num_valid = 0
        num_invalid = 0
        no_mapping = 0
        short_al = 0
        found_one = 0
        i = 0
        s1 = Set(self.p2cMapping.readToContig)
        s2 = Set(self.p2rMapping.readToContig)

        print len(self.p2cMapping.readToContig),\
        len(self.p2rMapping.readToContig), len(s1.intersection(s2))
        print len(self.p2cMapping.contigToRead)
        print len(self.c2rMapping.readToContig)
        print sum([len(c) for c in self.p2cMapping.readToContig.values()])
        print sum([len(c) for c in self.c2rMapping.readToContig.values()])
        print sum([len(c) for c in self.p2rMapping.readToContig.values()])
        for p in self.p2cMapping.readToContig:
            if i % 100 == 0:
                print i, num_valid, num_invalid, short_al,\
                no_mapping, found_one
            i += 1
            numfound = 0
            for pc in self.p2cMapping.readToContig[p]:
                p2c = PacbioAlignment(pc, self.p2cFileFormat)
                if not p2c.longAlignment(map_margin, p2c.seqLen, p2c.queryLen):
                    short_al += 1
                    continue
                found = 0
                if p2c.seqID in self.c2rMapping.readToContig and\
                p in self.p2rMapping.readToContig:
                    for cr in self.c2rMapping.readToContig[p2c.seqID]:
                        c2r = PacbioAlignment(cr, self.c2rFileFormat)
                        for pr in self.p2rMapping.readToContig[p]:
                            p2r = PacbioAlignment(pr, self.p2rFileFormat)
                            if self.valid_trio(p2c, c2r, p2r):
                                found = 1
                                break
                        if found == 1:
                            break
                    if found == 1:
                        num_valid += 1
                        numfound += 1
                    else:
                        num_invalid += 1
                else:
                    no_mapping += 1
            if numfound > 0:
                found_one += 1
        print num_valid, num_invalid, no_mapping

    def find_alignments(self, v, vStartMin, vStartMax, vDir):
        if v.Name() not in self.c2rMapping.readToContig:
            return []
        aligns = []
        for m in self.c2rMapping.readToContig[v.Name()]:
            a = PacbioAlignment(m, self.c2rFileFormat)
            if v.Strand() == '+':
                if vDir == '+':
                    if a.queryStrand != a.seqStrand:
                        continue
                    pos = a.seqStartPos - a.queryStartPos
                    if pos >= vStartMin and pos <= vStartMax:
                        aligns.append(a)
                else:
                    if a.queryStrand == a.seqStrand:
                        continue
                    pos = a.seqEndPos + a.queryStartPos
                    if pos >= vStartMin and pos <= vStartMax:
                        aligns.append(a)
            else:
                if vDir == '+':
                    if a.queryStrand == a.seqStrand:
                        continue
                    pos = a.seqStartPos - (a.queryLen - a.queryEndPos)
                    if pos >= vStartMin and pos <= vStartMax:
                        aligns.append(a)
                else:
                    if a.queryStrand != a.seqStrand:
                        continue
                    pos = a.seqEndPos + (a.queryLen - a.queryEndPos)
                    if pos >= vStartMin and pos <= vStartMax:
                        aligns.append(a)
        return aligns

    def path_from(self, path, pStartMin, pStartMax, pDir,
                  offset_margin=0.1):
        if len(path) == 0:
            return []
        start = path[0].v1
        e1 = path[0]
        a1s = self.find_alignments(start, pStartMin, pStartMax, pDir)
        if len(a1s) == 0:
            return None
        for a1 in a1s:
            offMin = (1 - offset_margin) * (start.length + e1.ovl)
            offMax = (1 + offset_margin) * (start.length + e1.ovl)
            if start.Strand() == '+':
                if pDir == '+':
                    sMin = a1.seqStartPos - a1.queryStartPos + offMin
                    sMax = a1.seqStartPos - a1.queryStartPos + offMax
                else:
                    sMax = a1.seqEndPos + a1.queryStartPos - offMin
                    sMin = a1.seqEndPos + a1.queryStartPos - offMax
            else:
                if pDir == '+':
                    sMin = a1.seqStartPos - (a1.queryLen - a1.queryEndPos) + \
                    offMin
                    sMax = a1.seqStartPos - (a1.queryLen - a1.queryEndPos) + \
                    offMax
                else:
                    sMax = a1.seqEndPos + (a1.queryLen - a1.queryEndPos) - \
                    offMin
                    sMin = a1.seqEndPos + (a1.queryLen - a1.queryEndPos) - \
                    offMax
            if len(path) > 1:
                p1 = self.path_from(path[1:], sMin, sMax, pDir, offset_margin)
                if p1 is None:
                    continue
                return [a1] + p1
            else:
                a2 = self.find_alignments(e1.v2, sMin, sMax, pDir)
                if len(a2) == 0:
                    continue
                return [a1, a2[0]]
        return None

    def valid_path(self, path, offset_margin=0.1):
        if len(path) == 0:
            return True
        if path[0].v1.Name() not in self.c2rMapping.readToContig:
            return None
        for m in self.c2rMapping.readToContig[path[0].v1.Name()]:
            a = PacbioAlignment(m, self.c2rFileFormat)
            if path[0].v1.Strand() == '+':
                if a.seqStrand == a.queryStrand:
                    pStart = a.seqStartPos - a.queryStartPos
                    p = self.path_from(path, pStart - 1, pStart + 1, '+',
                                       offset_margin)
                else:
                    pStart = a.seqEndPos + a.queryStartPos
                    p = self.path_from(path, pStart - 1, pStart + 1, '-',
                                       offset_margin)
            else:
                if a.seqStrand == a.queryStrand:
                    pStart = a.seqEndPos + (a.queryLen - a.queryEndPos)
                    p = self.path_from(path, pStart - 1, pStart + 1, '-',
                                       offset_margin)
                else:
                    pStart = a.seqStartPos - (a.queryLen - a.queryEndPos)
                    p = self.path_from(path, pStart - 1, pStart + 1, '+',
                                       offset_margin)
            if p is not None:
                return True
        return False

    def path_start_dir(self, path, offset_margin=0.1):
        if not self.valid_path(path, offset_margin):
            return None
        for m in self.c2rMapping.readToContig[path[0].v1.Name()]:
            a = PacbioAlignment(m, self.c2rFileFormat)
            if path[0].v1.Strand() == '+':
                if a.seqStrand == a.queryStrand:
                    pStart = a.seqStartPos - a.queryStartPos
                    p = self.path_from(path, pStart - 1, pStart + 1, '+',
                                       offset_margin)
                    pDir = '+'
                else:
                    pStart = a.seqEndPos + a.queryStartPos
                    p = self.path_from(path, pStart - 1, pStart + 1, '-',
                                       offset_margin)
                    pDir = '-'
            else:
                if a.seqStrand == a.queryStrand:
                    pStart = a.seqEndPos + (a.queryLen - a.queryEndPos)
                    p = self.path_from(path, pStart - 1, pStart + 1, '-',
                                       offset_margin)
                    pDir = '-'
                else:
                    pStart = a.seqStartPos - (a.queryLen - a.queryEndPos)
                    p = self.path_from(path, pStart - 1, pStart + 1, '+',
                                       offset_margin)
                    pDir = '+'
            if p is not None:
                return pStart, pDir
        return False
