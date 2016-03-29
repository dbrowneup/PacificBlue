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

#Author: Eric Fung
#Contact: kspham@eng.ucsd.edu

'''Simulated pacbio reads with C coverage
and read length follows a distribution H'''

import random
import bisect

with open('ecoli.fna') as f:
    header = f.readline()
    seq = f.read().replace('\n', '')
    seq_len = len(seq)

f.close()
print seq_len

coverage = 1000
maxlen = 20848
readlen = []

#create all the read length range
for i in xrange(0, maxlen, 500):
    value = []

    if i == 0:
        readlen.append(i + 500 / 2)
    else:
        readlen.append(i + 502 / 2)

#print readlen

prob = []

#assign prob
for i in xrange(len(readlen)):
    prob.append(.02439)

#print len(prob)
#print len(readlen)


#select the read lenght respect to the prob
def picklen(prob):
    breakpoints = [sum(prob[:x + 1]) for x in range(len(readlen))]
    #print bisect.bisect_left(breakpoints, random.uniform(0.0,breakpoints[-1]))
    return readlen[bisect.bisect_left(breakpoints,
                                      random.uniform(0.0, breakpoints[-1]))]

#crate all the read lenghts
with open('SimPacBioreads.txt', 'w') as f1:
    for i in xrange(coverage):
        start_loc = random.randint(0, seq_len)
        value = picklen(prob)

        if value % 2 != 1:
            end_loc = random.randint(value - 250, value + 250)
        else:
            end_loc = random.randint(value - 250, value + 249)

        if start_loc + 50 > seq_len:
            f1.write(seq[start_loc:seq_len] +
                     seq[0:(start_loc + end_loc) % seq_len] + '\n')

            print len(seq[start_loc:seq_len] +
                      seq[0:(start_loc + end_loc) % seq_len])
        else:
            f1.write(seq[start_loc:(start_loc + end_loc)] + '\n')

            print len(seq[start_loc:(start_loc + end_loc)])

f1.close()
