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

import sys


class AbstractAlignment:
    def __init__(self, aln, fileFormat="m4"):
        if fileFormat == "m5":
            self.qName = str(aln[0])
            self.qLength = int(aln[1])
            self.qStart = int(aln[2])
            self.qEnd = int(aln[3])
            self.qStrand = int(aln[4])
            self.tName = int(aln[5])
            self.tLength = int(aln[6])
            self.tStart = int(aln[7])
            self.tEnd = int(aln[8])
            self.tStrand = int(aln[9])
            self.score = int(aln[10])
            self.numMatch = int(aln[11])
            self.numMismatch = int(aln[12])
            self.numIns = int(aln[13])
            self.numDel = int(aln[14])
            self.mapQV = int(aln[15])
            self.qAlignedSeq = str(aln[16])
            self.matchPattern = str(aln[17])
            self.tAlignedSeq = str(aln[18])
        elif fileFormat == "m4":
            self.qName = str(aln[0])
            self.tName = int(aln[1])
            self.score = int(aln[2])
            self.percentSimilarity = float(aln[3])
            self.qStrand = int(aln[4])
            self.qStart = int(aln[5])
            self.qEnd = int(aln[6])
            self.qLength = int(aln[7])
            self.tStrand = int(aln[8])
            self.tStart = int(aln[9])
            self.tEnd = int(aln[10])
            self.tLength = int(aln[11])
            self.mapQV = int(aln[12])
        else:
            sys.exit("File format unrecognizable!")
    
    
    #ALL ORIGINAL CODE BELOW THIS STATEMENT - WILL NOT WORK RIGHT NOW (04/27/16 DB)
    #returns True if self contains align2 with a margin of epsilon
    def queryContains(self, align2, epsilon=1.0):
        if self.queryID != align2.queryID:
            return False
        if self.queryStartPos <= align2.queryStartPos and\
           self.queryEndPos >= align2.queryEndPos:
            if epsilon * (self.queryEndPos - self.queryStartPos) > \
                         (align2.queryEndPos - align2.queryStartPos):
                return True
        return False

    #returns True if self contains align2 with a margin of epsilon
    def seqContains(self, align2, epsilon=1.0):
        if self.seqID != align2.seqID:
            return False
        if self.seqStartPos <= align2.seqStartPos and \
           self.seqEndPos >= align2.seqEndPos:
            if epsilon * (self.seqEndPos - self.seqStartPos) > \
               (align2.seqEndPos - align2.seqStartPos):
                return True
        return False

    #returns true if alignment is sufficiently long or suffix prefix
    #in comparison to query/seq lengths
    def longAlignment(self, margin=0.5, seqLength=0, queryLength=0):
        #query contained in seq

        if seqLength == 0:
            seqLength = max(self.seqLen, self.seqEndPos)
        if queryLength == 0:
            queryLength = max(self.queryLen, self.queryEndPos)
        if margin * queryLength < self.queryEndPos - self.queryStartPos:
            return True
        #seq contained in query
        if margin * seqLength < self.seqEndPos - self.seqStartPos:
            return True
        #prefix of query
        if self.queryStartPos < (1 - margin) * \
                                (self.queryEndPos - self.queryStartPos):
            #suffix of sequence
            if self.seqStrand == self.queryStrand:
                if (seqLength - self.seqEndPos) < \
                   (1 - margin) * (self.seqEndPos - self.seqStartPos):
                    return True
            #prefix of sequence
            else:
                if self.seqStartPos < \
                        (1 - margin) * (self.seqEndPos - self.seqStartPos):
                    return True
        #suffix of query
        if (queryLength - self.queryEndPos) < \
                (1 - margin) * (self.queryEndPos - self.queryStartPos):
            #suffix of sequence
            if self.seqStrand != self.queryStrand:
                if (seqLength - self.seqEndPos) < \
                        (1 - margin) * (self.seqEndPos - self.seqStartPos):
                    return True
            #prefix of sequence
            else:
                if self.seqStartPos < \
                        (1 - margin) * (self.seqEndPos - self.seqStartPos):
                    return True

        return False

    def ConcordantDirection(self, align2, v1, v2):
        a1 = self
        a2 = align2
        if (a1.queryID != a2.queryID):
            sys.exit("Inconsistent mapping and vertices")
        if a1.seqID != v1.Name() or a2.seqID != v2.Name():
            sys.exit("Inconsistent mapping and vertices")

        if (a1.seqStrand == a1.queryStrand):
            strand1 = '+'
        else:
            strand1 = '-'

        if (a2.seqStrand == a2.queryStrand):
            strand2 = '+'
        else:
            strand2 = '-'

        #outgoing for one vertex and incoming for another
        if (strand1 == strand2) and v1.Strand() == v2.Strand():
            return True
        elif (strand1 != strand2) and v1.Strand() != v2.Strand():
            return True
        else:
            return False

    def relative_offset(self, align2, v1, v2):
        a1 = self
        a2 = align2
        if a1.seqID != v1.Name() or a2.seqID != v2.Name():
            return False
        if a1.queryID != a2.queryID:
            return False
        if not a1.ConcordantDirection(align2, v1, v2):
            return False

        toggleP1 = False
        toggleP2 = False

        if v1.Strand() == '+':
            s1 = a1.seqStartPos
        else:
            s1 = v1.length - a1.seqEndPos
            toggleP1 = not toggleP1
        if v2.Strand() == '+':
            s2 = a2.seqStartPos
        else:
            s2 = v2.length - a2.seqEndPos
            toggleP2 = not toggleP2

        if a1.queryStrand == '-':
            toggleP1 = not toggleP1
        if a2.queryStrand == '-':
            toggleP2 = not toggleP2

        if a1.seqStrand == '-':
            toggleP1 = not toggleP1
        if a2.seqStrand == '-':
            toggleP2 = not toggleP2
        #outgoing edge
        if not toggleP1:
            p1 = a1.queryStartPos
        else:
            p1 = a1.queryLen - a1.queryEndPos
        #incoming edge
        if not toggleP2:
            p2 = a2.queryStartPos
        else:
            p2 = a2.queryLen - a2.queryEndPos

        return ((p2 - s2) - (p1 - s1))

    def ConcordantEdge(self, align2, edge, offset_margin=0.1):
        a1 = self
        a2 = align2
        if a1.seqID != edge.v1.Name() or a2.seqID != edge.v2.Name():
            return False
        if a1.queryID != a2.queryID:
            return False
        if not a1.ConcordantDirection(align2, edge.v1, edge.v2):
            return False
        l1 = edge.v1.length
        e_offset = l1 + edge.ovl
        r_offset = self.relative_offset(align2, edge.v1, edge.v2)

        if (abs(float(r_offset - e_offset)) / e_offset < offset_margin):
            return True
        return False

    def spansEdge(self, align2, edge, offset_margin=0.1, span_margin=100):
        a1 = self
        a2 = align2
        if not a1.ConcordantEdge(a2, edge, offset_margin):
            return False
        if a1.seqStrand == a1.queryStrand:
            if a1.seqStartPos >= edge.v1.length + (edge.ovl - span_margin) - 1:
                return False
        else:
            if a1.seqEndPos <= -1 * (edge.ovl - span_margin) + 1:
                return False
        if a2.seqStrand == a2.queryStrand:
            if a2.seqEndPos <= -1 * (edge.ovl - span_margin) + 1:
                return False
        else:
            if a2.seqStartPos >= edge.v2.length + (edge.ovl - span_margin) - 1:
                return False
        return True

    def align_length(self):
        return self.queryEndPos - self.queryStartPos
