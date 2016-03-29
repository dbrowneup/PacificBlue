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
from path.path import path


class PacbioMapping:

    def __init__(self, fileName, fileFormat="m4", map_margin=0.5):
        self.readToContig = {}
        self.contigToRead = {}
        self.alignments = []
        self.numAlignments = 0
        self.fileFormat = fileFormat

        if fileFormat in ['m4', 'blastn', 'm4algae']:
            fpath = path(fileName)
            f = fpath.open()
            for line in f:
                if len(line) <= 1 or line[0] == "#" or \
                    line == "qname tname score pctsimilarity qstrand " + \
                            "qstart qend qseqlength tstrand tstart tend " + \
                            "tseqlength mapqv ncells clusterScore probscore" \
                            + " numSigClusters\n" or \
                    line == "qname tname qstrand tstrand score " + \
                            "pctsimilarity tstart tend tlength qstart qend" + \
                            " qlength ncells\n":
                    continue
                self.alignments.append(line.strip())
                self.numAlignments += 1
                pa = PacbioAlignment(line.strip(), fileFormat)
                if not pa.longAlignment(map_margin):
                    continue
#                 if fileFormat == 'blastn':
#                     ll = line.strip().split()
#                     if float(ll[2]) < 98:
#                         continue
                sID = pa.seqID
                qID = pa.queryID
                if sID in self.contigToRead:
                    self.contigToRead[sID]. \
                        append(self.alignments[self.numAlignments - 1])
                else:
                    self.contigToRead[sID] = \
                        [self.alignments[self.numAlignments - 1]]

                if qID in self.readToContig:
                    self.readToContig[qID]. \
                        append(self.alignments[self.numAlignments - 1])
                else:
                    self.readToContig[qID] = \
                        [self.alignments[self.numAlignments - 1]]
            f.close()
        #sort mappings in decreasing order of contig length
        for qID in self.readToContig:
            self.readToContig[qID].sort(cmp=lambda x, y:
                                        PacbioAlignment(y, fileFormat).seqLen -
                                        PacbioAlignment(x, fileFormat).seqLen)

    def positions_on_seq(self, qID, qStrand):
        if qID not in self.readToContig:
            return []
        positions = []
        for m in self.readToContig[qID]:
            a = PacbioAlignment(m, self.fileFormat)
            if qStrand == '+':
                if a.queryStrand == a.seqStrand:
                    positions.append(a.seqStartPos - a.queryStartPos)
                else:
                    positions.append(a.seqEndPos + a.queryStartPos)
            else:
                if a.queryStrand == a.seqStrand:
                    positions.append(a.seqEndPos +
                                     (a.queryLen - a.queryEndPos))
                else:
                    positions.append(a.seqStartPos -
                                     (a.queryLen - a.queryEndPos))
        return positions
#m = PacbioMapping("data/sorted_best30.fasta.m4")
#print m.contigToRead["228045"]
