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

import experimental
from PacbioAlignment import PacbioAlignment
from validation.mapping_validation import MappingValidation
from PacbioMapping import PacbioMapping
#
mv = MappingValidation(experimental.contig_refmap_file, experimental.contig_pacbiomap_file, experimental.pacbio_refmap_file)
m = PacbioMapping(experimental.contig_pacbiomap_file)
#
#count = {}
#    
#x = {x:len(mv.c2rMapping.contigToRead[x]) for x in mv.c2rMapping.contigToRead.keys()}
#

#y = {x:PacbioAlignment(mv.c2rMapping.readToContig[x][0],"blastn").queryStartPos for x in mv.c2rMapping.readToContig.keys()}

#
#xkeys = y.keys()

#xkeys.sort(lambda b,q: y[b]-y[q])
#p = [(d, x[d], y[d]) for d in xkeys]
#for d in p:
#    print d[0], d[1], d[2]

def count_read_map(contig_id):
    temp = []
    mapp = []
    if str(abs(contig_id)) in m.contigToRead.keys():
        for i in xrange(len(m.contigToRead[str(abs(contig_id))])):
            pa = PacbioAlignment(m.contigToRead[str(abs(contig_id))][i],'m4')
            paID = pa.queryID
            if pa not in temp:            
                temp.append(pa)
                for x in m.readToContig[paID] :
                    mapp.append(x)
                           
        return temp 
for j in  (set(count_read_map(734)) & set(count_read_map(714))):
    print j


def get_mapping_from_read(reads):
    mapp = []
    for r in reads:
        if r in m.readToContig:
            for k in m.readToContig[r]:
                mapp.append(k)
    return mapp


reads = [40797, 63328,19978,34685,60665,57663,73914,35859,18131,65893,45985,28350,31858,42635,1552,42852]

#for j in get_mapping_from_read(reads):
#    print j

class Find_reads(object):
    def __init__(self, gap_id):
        self.gap_id = gap_id
        self.count = 0
        self.reads = []

    def add_read(self, read):        
        if read not in self.reads:
            self.reads.append(read)
            self.count += 1
#            print self.count

    def get_gap_id(self):
        return self.gap_id

    def get_count(self):
        return self.count

    def __repr__(self):
        return "%s" % ( self.get_count())

def find_gaps(start, end):
    m = PacbioMapping(experimental.pacbio_refmap_file, 'm4', .7)
    mapping = []
    count = {} 

    for x in m.readToContig:
        for y in m.readToContig[x]:
            pa = PacbioAlignment(y, 'm4')
            pasta = pa.seqStartPos
            paend = pa.seqEndPos

        if pasta < start and paend > end:
#            if 2 not in count:
#                count[2] = 0
#            count[2] += 1
           if x not in mapping: 
                mapping.append(x)

    return mapping

def run():
    len_683 = 1702
    len_670 = 2463
    len_746 = 2550 
    
    gaps = [[225327 - len_683, 228688],
            [ 732362 - len_746,  735306  ],     
            [ 1430588,  1432760 ],
            [  1633434 - len_670, 1632578],    
            [ 2726495,  2727728 + len_683],
            [ 3421800,  3425228 + len_683],
            [ 3617367,  3619916],   
            [ 3763506- len_746,  3763594],
            [ 3941387 - len_683,  3944643],
            [ 4035109 - len_683,  4038795],
            [ 4166237 - len_683,  4169926],
            [ 4207726 - len_683,  4209938]]
    
    
    a = []
    for i,q in enumerate(gaps):
    #    print i+1
        for w in find_gaps(q[0], q[1]):
            a.append(w)
     #   print " \n"
    for j, x in enumerate(get_mapping_from_read(a)):
        
        print x
        
#run()    

#        if pasta < 225327 - 1702 and paend > 228688:
#            a = Find_reads(1) 
##            print a
#            
#            if a not in mapping.values():
#                mapping[1] = a
#            mapping[1].add_read(y) 
##            if y not in mapping.values():
##                mapping.values().append(y)
##
##            if 1 not in count:
##                count[1] = 0
##            count[1] += 1
#        elif pasta < 732362 - 2550 and paend > 735306:
#            b = Find_reads(2) 
#            if b not in mapping.values():
#                mapping[2] = (b)
#            mapping[2].add_read(y) 
##
##            if y not in mapping.values():
##                mapping.values().append(y)  
##         
##            if 2 not in count:
##                count[2] = 0
##            count[2] += 1
#
#        elif pasta < 1430588 and paend > 1432760 - 2463:
#
#            c = Find_reads(3) 
#            if c not in mapping.values():
#                mapping[3] = (c)
#            mapping[3].add_read(y) 
##
#
#
##            if y not in mapping.values():
##                mapping.values().append(y)
##          
##            if 3 not in count:
##                count[3] = 0
##            count[3] += 1
#
#        elif pasta < 1632578 and paend > 1634343 - 2463:
#            d = Find_reads(4) 
#            if d not in mapping.values():
#                mapping[4] = d
#            mapping[4].add_read(y) 
##
##            if y not in mapping.values():
##                mapping.values().append(y)
##           
##            if 4 not in count:
##                count[4] = 0
##            count[4] += 1
#
#        elif pasta < 2726495 and paend > 2727728:
#            d = Find_reads(5) 
#            if d not in mapping.values():
#                mapping[5] = d
#            mapping[5].add_read(y) 
##       
##            if y not in mapping.values():
##                mapping.values().append(y) 
##
##            if 5 not in count:
##                count[5] = 0
##            count[5] += 1
#
#        elif pasta < 3421800 and paend > 3425228:
#            d = Find_reads(6) 
#            if d not in mapping.values():
#                mapping[6] = d
#            mapping[6].add_read(y) 
##
##            if y not in mapping.values():
##                mapping.values().append(y)
##
##            if 6 not in count:
##                count[6] = 0
##            count[6] += 1
#
#        elif pasta < 3619916 - 2550 and paend > 3617367:
#            d = Find_reads(7) 
#            if d not in mapping.values():
#                mapping[7] = (d)
#            mapping[7].add_read(y) 
##       
##            if y not in mapping.values():
##                mapping.values().append(y) 
##            if 7 not in count:
##                count[7] = 0
##            count[7] += 1
#
#
#        elif pasta < 3763506 - 2550 and paend > 3763594:
#            d = Find_reads(8) 
#            if d not in mapping.values():
#                mapping[8] = (d)
#            mapping[8].add_read(y) 
##       
##            if y not in mapping.values():
##                mapping.values().append(y)
##
##            if 8 not in count:
##                count[8] = 0
##            count[8] += 1
#
#        elif pasta < 3941387 - 1702 and paend > 3944643:
#            d = Find_reads(9) 
#            if d not in mapping.values():
#                mapping[9] = d
#            mapping[9].add_read(y) 
##
##            if y not in mapping.values():
##                mapping.values().append(y)
##
##            if 9 not in count:
##                count[9] = 0
##            count[9] += 1
#       
#        elif pasta < 4035109 - 1702 and paend > 4038795:
#            d = Find_reads(10) 
#            if d not in mapping.values():
#                mapping[10] = (d)
#            mapping[10].add_read(y) 
##
##            if y not in mapping.values():
##                mapping.values().append(y) 
##
##            if 10 not in count:
##                count[10] = 0
##            count[10] += 1
#       
#        elif pasta < 4166237  -1702 and paend > 4169926:
#            d = Find_reads(11) 
#            if d not in mapping.values():
#                mapping[11] = d
#            mapping[11].add_read(y) 
##       
##            if y not in mapping.values():
##                mapping.values().append(y)
##
##            if 11 not in count:
##                count[11] = 0
##            count[11] += 1
##       
#        elif pasta < 4207726 - 1702 and paend > 4209937:
#            d = Find_reads(12) 
#
#            if d not in mapping.values():
#                mapping[12] = d
#            mapping[12].add_read(y) 
##
##            if y not in mapping.values():
##                mapping.values().append(y)
##
##            if 12 not in count:
##                count[12] = 0
##            count[12] += 1
##       
#print (mapping)
#mapping.values().sort(lambda x, y: PacbioAlignment(x).seqStartPos - PacbioAlignment(y).seqStartPos)


#for x in mapping.values():
#    print x
#    for i in xrange(len(m.contigToRead[j])):
#        pa = PacbioAlignment(m.contigToRead[j][i],'m4')
#        pasta = pa.queryStartPos
#        print pasta

 
#for x in m.readToContig:
#pa = PacbioAlignment(m, 'm4')
#print pa.queryStartPos
#
