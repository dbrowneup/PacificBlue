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

from AbstractAlignment import AbstractAlignment


class PacbioAlignment(AbstractAlignment):

    def __init__(self, line, fileFormat='m4'):
        def Strand(ss):
            if ss == '0':
                return '+'
            return '-'

        ll = line.strip().split()
        if fileFormat == 'm4':
            if Strand(ll[3]) == '-':
                sStart = int(ll[8]) - int(ll[7])
                sEnd = int(ll[8]) - int(ll[6])
            else:
                sStart = int(ll[6])
                sEnd = int(ll[7])
            if Strand(ll[2]) == '-':
                qStart = int(ll[11]) - int(ll[10])
                qEnd = int(ll[11]) - int(ll[9])
            else:
                qStart = int(ll[9])
                qEnd = int(ll[10])
            AbstractAlignment.__init__(self, ll[1], Strand(ll[3]),
                                       sStart, sEnd, ll[0],
                                       Strand(ll[2]), qStart, qEnd,
                                       ll[8], ll[11])
        elif fileFormat == 'm4algae':
            if Strand(ll[8]) == '-':
                sStart = int(ll[11]) - int(ll[10])
                sEnd = int(ll[11]) - int(ll[9])
            else:
                sStart = int(ll[9])
                sEnd = int(ll[10])
            if Strand(ll[4]) == '-':
                qStart = int(ll[7]) - int(ll[6])
                qEnd = int(ll[7]) - int(ll[5])
            else:
                qStart = int(ll[5])
                qEnd = int(ll[6])
            AbstractAlignment.__init__(self, ll[1], Strand(ll[8]), sStart,
                                       sEnd, ll[0], Strand(ll[4]),
                                       qStart, qEnd, ll[11], ll[7])
        elif fileFormat == 'blastn':
            sStrand = '+'
            if len(ll) >= 14:
                if int(ll[9]) - int(ll[8]) < 0:
                    sStrand = '-'
                    AbstractAlignment.__init__(self, ll[1], sStrand, ll[9],
                                               ll[8], ll[0], '+', ll[6], ll[7],
                                               ll[13], ll[12])
                else:
                    sStrand = '+'
                    AbstractAlignment.__init__(self, ll[1], sStrand, ll[8],
                                               ll[9], ll[0], '+', ll[6], ll[7],
                                               ll[13], ll[12])
            else:
                if int(ll[9]) - int(ll[8]) < 0:
                    sStrand = '-'
                    AbstractAlignment.__init__(self, ll[1], sStrand, ll[9],
                                               ll[8], ll[0], '+', ll[6], ll[7])
                else:
                    sStrand = '+'
                    AbstractAlignment.__init__(self, ll[1], sStrand, ll[8],
                                               ll[9], ll[0], '+', ll[6], ll[7])

        else:
            ll = line.strip().split()
            AbstractAlignment.__init__(self, ll[0], ll[1], ll[2],
                                       ll[3], ll[4], ll[4], ll[5], ll[6])
