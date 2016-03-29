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

import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
matplotlib.rcParams['interactive'] = "true"


class ReadLengthStatistics:

    def __init__(self, faFile):
        self.readLen = {}
        self.faFile = faFile
        f = open(self.faFile)
        cl = 0
        name = ''
        for line in f:
            if len(line.strip()) == 0 or line[0] == '#':
                continue
            if line[0] == '>':
                if cl != 0:
                    self.readLen[name] = cl
                    cl = 0
                name = line[1:].strip().split()[0]
                continue
            cl += len(line.strip())
        f.close()
        if cl != 0:
            self.readLen[name] = cl
            cl = 0

    def LengthHistogram(self, pngFile='test.png',
                        includeCondition=lambda x, y: True):
        lens = []
        nr = 0
        for r in self.readLen:
            nr += 1
#             if nr % 10000 == 0:
#                 print nr
            if includeCondition(r, self.readLen[r]):
                lens.append(self.readLen[r])
        lens.sort()

        def getNxx(lens, eps):
            ss = sum(lens)
            st = 0
            for l in lens:
                st += l
                if st >= ss * (1 - eps):
                    return l
            return lens[-1]

        print "#=============================================================="
        print "#", self.faFile
        print "#Num seq =", len(lens)
        print "#Max seq length =", max(lens)
        print "#Min seq length =", min(lens)
        print "#Average seq length =", \
            (lambda x: 0 if len(x) == 0 else sum(x) / len(x))(lens)
        print "#Median seq length =", \
            (lambda x: 0 if len(x) == 0 else x[len(x) / 2])(lens)
        print "#N50 seq length =", getNxx(lens, 0.5)
        print "#N90 seq length =", getNxx(lens, 0.9)
        print "#=============================================================="

        maxl = max(lens)
        minl = min(lens)
        bins = {2 ** i: 0 for i in range(int(math.log(maxl, 2) + 1))}
        for cl in lens:
            bins[2 ** int(math.log(cl, 2))] += 1
        binkeys = bins.keys()
        binkeys.sort()
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        ind = np.arange(int(math.log(maxl, 2)) + 1 - int(math.log(minl, 2)))
        width = 0.7
        p1 = plt.bar(ind, [bins[bk] for bk in \
                binkeys[int(math.log(minl,2)):]], width, color='r')
        plt.ylabel('Number of queries')
        plt.title('Length distribution of ' + self.faFile)
        plt.xlabel('Length of query')
        plt.xticks(ind + width / 2, binkeys[int(math.log(minl, 2)):])
        #plt.legend((p1[0]), ('Pacbio Read Lengths')))
        plt.savefig(pngFile)
        plt.show()


#(ReadLengthStatistics('data/spalgae/spalgae-velvet.fa')).\
#    LengthHistogram(pngFile='data/spalgae/spalgae-velvet.png')
