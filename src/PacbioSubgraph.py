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

#Re-written by Dan Browne on 04/27/16

from illumina_graph import Graph
from illumina_path import IlluminaPath
from utils import toggleStrand
import numpy as np
import itertools

class PacbioSubgraph(Graph):

    def __init__(self, pacbio_id, illumina_graph, pacbio_mapping):
        Graph.__init__(self)
        self.pacbio_id = pacbio_id
        self.FivePOvhg = []
        self.ThrePOvhg = []
        self.Connects = []
        #Create array of read aligments by read length
        self.pacbio_mapping = pacbio_mapping
        self.readArray = self.read_stack()
        qStart = 0
        qEnd = 0
        #Mark coordinates of alignments on read
        for n in range(len(self.pacbio_mapping.readToContig[pacbio_id])):
            qStart = self.pacbio_mapping.readToContig[pacbio_id][n].qStart
            qEnd = self.pacbio_mapping.readToContig[pacbio_id][n].qEnd
            self.mark_coords(n, qStart, qEnd)
        #Calculate per-base coverage of read
        self.covArray = self.readArray.sum(axis=0)
        print "Total coverage of read:", self.covArray.sum()
        #Filter out alignments in repetitive regions
        align_coords = set([])
        repeat_coords = set(self.find_repeats())
        for n in self.pacbio_mapping.readToContig[pacbio_id]:
            align_coords = set(range(n.qStart, n.qEnd))
            self.filter_repeats(align_coords, repeat_coords, n)
        #Find 5' and 3' overhangs
        for n in self.pacbio_mapping.readToContig[pacbio_id]:
            self.find_overhangs(n)
        #Connect 5' and 3' overhangs
        possible_connects = list(itertools.product(self.ThrePOvhg, self.FivePOvhg))
        self.Connects = map(self.test_connect, possible_connects)
        
    
    def read_stack(self):
        stack_size = len(self.pacbio_mapping.readToContig[self.pacbio_id])
        read_size = self.pacbio_mapping.readToContig[self.pacbio_id][0].qLength
        a = np.zeros((stack_size, read_size))
        print "Shape of array:", np.shape(a)
        return a
    
    def mark_coords(self, n, qStart, qEnd):
        for i in range(qStart, qEnd):
            self.readArray[n][i] = 1
    
    def find_repeats(self, cutoff=3):
        repeat_coords = []
        for i in range(len(self.covArray)):
            if self.covArray[i] > cutoff:
                repeat_coords.append(i)
        print str(len(repeat_coords))+" repetitive bps"
        return repeat_coords

    def filter_repeats(self, align_coords, repeat_coords, n, fraction=0.5):
        bp_in_repeat = align_coords.intersection(repeat_coords)
        alignment = self.pacbio_mapping.readToContig[self.pacbio_id]
        if len(bp_in_repeat) >= fraction*len(align_coords):
            alignment = alignment.remove(n)

    def find_overhangs(self, a):
        if a.qStart > a.tStart and (a.qLength - a.qEnd) < (a.tLength - a.tEnd):
            print "5' overhang"
            self.FivePOvhg.append(a)
        elif a.qStart < a.tStart and (a.qLength - a.qEnd) > (a.tLength - a.tEnd):
            print "3' overhang"
            self.ThrePOvhg.append(a)
        elif a.qStart > a.tStart and (a.qLength - a.qEnd) > (a.tLength - a.tEnd):
            print "Double overhang"
        elif a.qStart <= a.tStart and (a.qLength - a.qEnd) <= (a.tLength - a.tEnd):
            print "No overhang"
        else:
            print "Unknown overhang!"
    
    def test_connect(self, p):
        tp, fp = p
        if tp.qEnd < fp.qStart:
            if (tp.tLength - tp.tEnd + fp.tStart) < (fp.qStart - tp.qEnd):
                print "Read spans gap between scaffolds"
            else:
                print "Gap between alignments, possible overlap of scaffold ends"
        elif fp.qEnd >= tp.qEnd >= fp.qStart >= tp.qStart:
            if fp.tEnd < (fp.qEnd - tp.qStart) > (fp.tLength - fp.tStart):
                print "Read spans overlap of scaffold ends"
            else:
                print "Scaffolds extend beyond alignment overlap"
        elif tp.qStart > fp.qEnd:
            print "Gap between inverted alignments, unlikely scaffold overlap"
        else:
            print "Strange alignment behavior"
            print "Alignment coordinates (qStart, qEnd):"
            print "A1 (%d, %d)" % (tp.qStart, tp.qEnd)
            print "A2 (%d, %d)" % (fp.qStart, fp.qEnd)
