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

#Re-written by Dan Browne on 04/05/16

from AbstractAlignment import AbstractAlignment
from datetime import datetime
import matplotlib
import matplotlib.pyplot as plt

class PacbioMapping:

    def __init__(self, fileName, fileFormat="m4"):
        print "Beginning to load PacBio mapping:", str(datetime.now())
        self.readToContig = {}
        self.fileFormat = fileFormat
        self.alignments = open(fileName, "rU").read().split('\n') # reads alignments
        self.alignments = self.alignments[1:-1] if "score" in self.alignments[0] else self.alignments[0:-1] # discards empty last row and first row if header
        self.alignments = [x.split(' ') for x in self.alignments] # splits alignment into list
        self.readToContig = {x[0]: [] for x in self.alignments if x[0] not in self.readToContig} # parses qname into dictionary
        #parse alignments into dictionaries
        for align in self.alignments:
            self.readToContig[align[0]].append(AbstractAlignment(align, fileFormat))
        print "Number of mapping reads:", len(self.readToContig)
        print "Number of alignments:", len(self.alignments)
        #sort readToContig mappings in decreasing order of contig length
        for qID in self.readToContig:
            self.readToContig[qID] = sorted(self.readToContig[qID], key=lambda x: x.tLength, reverse=True)
        print "Finished loading PacBio mapping:", str(datetime.now())
    def filter_reads(self, length_fraction=0.7):
        #remove all reads with only 1 alignment
        #remove reads where 1st alignment covers >= 70% (length_fraction) of read length
        print "Beginning read filtration:", str(datetime.now())
        filtered_reads = []
        for k,v in self.readToContig.items():
            map_length = v[0].qEnd - v[0].qStart
            if len(v) == 1:
                filtered_reads.append(k)
                del self.readToContig[k]
            elif map_length >= v[0].qLength * float(length_fraction):
                filtered_reads.append(k)
                del self.readToContig[k]
        self.alignments = [x for x in self.alignments if x[0] not in filtered_reads]
        print "Number of mapping reads filtered out:", len(filtered_reads)
        print "Number of mapping reads remaining:", len(self.readToContig)
        print "Finished read filtration:", str(datetime.now())
    def read_mapping_frequency(self, title, bins=30):
        x_axis = [len(x) for x in self.readToContig.values()]
        plt.hist(x_axis, bins)
        plt.xlabel("Alignments per read")
        plt.ylabel("Number of reads")
        plt.savefig(title+".png")
        plt.clf()

