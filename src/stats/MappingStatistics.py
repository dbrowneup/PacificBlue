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

#import numpy as np
#import matplotlib.pyplot as plt
#from sets import Set
#import math
import matplotlib
matplotlib.rcParams['interactive'] = "true"

#from PacbioAlignment import PacbioAlignment
from illumina_graph import Graph
from PacbioMapping import PacbioMapping
from PacbioSubgraph import PacbioSubgraph

from ReadLengthStatistics import ReadLengthStatistics


class MappingStatistics:

    def __init__(self, queryFa, refFa, mappingFileName,
                 fileFormat='m4', contigDot=''):
        print "#loading pacbio mapping"
        self.mapping = PacbioMapping(mappingFileName, fileFormat)
        print "#pacbio mapping loaded"
        self.queryFa = queryFa
        self.refFa = refFa
        self.queryLenStats = ReadLengthStatistics(self.queryFa)
        self.refLenStats = ReadLengthStatistics(self.refFa)
        self.contigDot = contigDot
        print "#mapping statistics initialized"

    def QueryNumMapping(self, map_margin=0.5, offset_margin=0.1):
        #nn = [len(i) for i in [j for j in self.mapping.readToContig.values() \
            #  if PacbioAlignment(j).longAlignment()]]

#        nn = []
#        rr = 0
#        rt = 0
#        ll = {}
#        for j in self.mapping.readToContig:
#            n = 0
#            ll[j] = 0
#            for i in self.mapping.readToContig[j]:
#                k = PacbioAlignment(i)
#                self.refLenStats.readLen[k.seqID]
#                if k.queryID not in self.queryLenStats.readLen:
#                    rr += 1
#                    continue
#                else:
#                    rt += 1
#                if k.longAlignment(self.refLenStats.readLen[k.seqID], \
#                   self.queryLenStats.readLen[k.queryID], margin):
#                    n += 1
#                    ll[j] += 1
#            nn.append(n)
        #print rr, rt
        #bins = [0 for i in range(max(nn)+1)]
        #for n in nn:
        #    bins[n] += 1
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ind = np.arange(max(nn)+1)
        #width = 0.7
        #p1 = plt.bar(ind, bins, width, color='r')
        #plt.ylabel('Number of queries')
        #plt.title('Distribution of reported mappings of \n' + \
            #self.queryFa+ ' to ' + self.refFa +\
            #'\n with match length at least '+str(100*margin)+' percent')
        #plt.xlabel('Number of mappings query')
        #plt.xticks(ind+width/2, range(max(nn)+1))
        ##plt.legend((p1[0]), ('Pacbio Read Lengths')))
        #plt.savefig("data/spalgae/spalgae_pacbio_illumina"+\
                        #.NumMapping_len50.png")
        #self.queryLenStats.LengthHistogram(pngFile=\
        #'data/spalgae/spalgae_pacbio_illumina_mapping50_length.png', \
        #includeCondition = lambda x, y: True if x in \
        #self.mapping.readToContig and len(self.mapping.readToContig[x])>0 \
        #and ll[x] > 0 else False)
        #self.queryLenStats.LengthHistogram(pngFile=\
        #'data/spalgae/spalgae_pacbio_illumina_2mapping50_length.png', \
        #includeCondition = lambda x, y: True if x in \
        #self.mapping.readToContig and len(self.mapping.readToContig[x])>1 \
        #and ll[x] > 1 else False)
        #self.queryLenStats.LengthHistogram(pngFile=\
        #'data/spalgae/spalgae_pacbio_illumina_3mapping50_length.png',\
        #includeCondition = lambda x, y: True if x in\
        #self.mapping.readToContig and len(self.mapping.readToContig[x])>2\
        #and ll[x] > 2 else False)

        if self.contigDot != '':
            contigGraph = Graph()
            contigGraph.load(self.contigDot, self.refFa)
            print "#Graph loaded"

        def graphHasEdge(x, y):
            return True
            if ((x not in self.mapping.readToContig or
                 len(self.mapping.readToContig[x]) < 1)):  # or ll[x] <= 1:
                    return False
            sg = PacbioSubgraph(x, contigGraph, self.mapping,
                                map_margin, True, offset_margin)
            lvs = len(sg.vs)
            if lvs > 2:
                return len(sg.es) > 0
            return len(sg.es) > 0

        self.queryLenStats.LengthHistogram(pngFile='abc.png',
            #pngFile='../data/SAureusTW20/' +
            #'SAureusTW20_pacbio_miseqcontig_mappingedge.png',
            includeCondition=graphHasEdge)
