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

from stats.MappingStatistics import MappingStatistics
from validation.mapping_validation import MappingValidation
from long_contig_graph import LongContigGraph

from illumina_graph import Graph
from PacbioMapping import PacbioMapping
from PacbioAlignment import PacbioAlignment
from sets import Set


dataName = 'ecoli0'
pbMappingFormat = 'm4'
c2rMappingFormat = 'blastn'
p2cMappingFormat = 'm4'
p2rMappingFormat = 'm4'

dataDir = '../data/'
prefix = dataDir + dataName + '/' + dataName
contigFa = prefix + '-contigs.fa'
contigDot = prefix + '-contigs.dot'
pbMapping = prefix + '_pacbio_contigs_mapping.fasta.m4'
outFa = prefix + '_long_contig_extended.fa'
outDot = prefix + '_long_contig_extended.dot'

c2rMapping = prefix + '_contigs_ref_mapping.blastn'
p2cMapping = prefix + '_pacbio_contigs_mapping.fasta.m4'
p2rMapping = prefix + '_pacbio_ref_mapping.fasta.m4'

min_vertex_length = 2048


#
#
# m = MappingStatistics('../data/spalgae/spalgae_pacbio.fa',
#                       '../data/spalgae/spalgae-contigs.fasta',
#                       '../data/spalgae/sorted_best30.fasta.m4',
#                       'm4algae', '../data/spalgae/spalgae-contigs.dot'
#                       ).QueryNumMapping()
#
#
#

# m = MappingStatistics('../data/ecoli/ecoli_pacbio.fa',
#                      '../data/ecoli/ecoli-contigs.fa',
#                      '../data/ecoli/ecoli_pacbio_contigs_mapping.fasta.m4',
#                      'm4', '../data/ecoli/ecoli-contigs.dot'
#                      ).QueryNumMapping()
# quit()

# m = MappingStatistics('../data/SAureusTW20/SAureusTW20_pacbio.fa',
#                     '../data/SAureusTW20/SAureusTW20_miseq-contigs.fa',
#                     '../data/SAureusTW20/SAureusTW20_pacbio_miseqcontig_mapping.fasta.m4',
#                     'm4', '../data/SAureusTW20/SAureusTW20_miseq-contigs.dot'
#                      )
# m.QueryNumMapping()
# m.refLenStats.LengthHistogram(pngFile='../data/SAureusTW20/' +
#                               'SAureusTW20_miseqcontig_length.png')
#
#
#m = MappingStatistics('../data/test/path_pacbio.fa',
#                      '../data/test/path_contigs.fa',
#                      '../data/test/path_pacbio_contig.fasta.m4',
#                      'm4', '../data/test/path_contigs.dot'
#                      ).QueryNumMapping()
#
#r2c = m.mapping.readFilteredMappings(m.queryLenStats.readLen,
#                                     m.refLenStats.readLen)
#g = open("data/ecoli_contigs_ref_filtered.blastn", "w")
#for r in r2c:
#    for a in r2c[r]:
#        g.write(a + '\n')
#g.close()
#
mv = MappingValidation(c2rMapping, p2cMapping, p2rMapping, c2rMappingFormat,
                       p2cMappingFormat, p2rMappingFormat)

# mv = MappingValidation('../data/SAureusTW20/SAureusTW20_miseqcontigs_ref_mapping.blastn',
#                        '../data/SAureusTW20/SAureusTW20_pacbio_miseqcontig_mapping.fasta.m4',
#                        '../data/SAureusTW20/SAureusTW20_pacbio_ref_mapping.fasta.m4')
# print len(mv.c2rMapping.readToContig)
# print len([p for p in mv.p2cMapping.readToContig if\
#            len(mv.p2cMapping.readToContig[p]) > 1])
# mv.trio_stats()

ig = Graph()
ig.load(contigDot, contigFa)
print 'Graph loaded', len(ig.vs), len(ig.es)

# f = open(dataDir+dataName+'/'+dataName+'_length_coverage.stats', "w")
# for vid in ig.vs:
#     if vid < 0:
#         continue
#     v = ig.vs[vid]
# #     if str(vid) not in mv.c2rMapping.readToContig:
# #         print v.length, v.cov/(v.length - 63), 0
# #     else:
# #         print v.length, v.cov/(v.length - 63), len(mv.c2rMapping.readToContig[str(vid)])
#     f.write(str(vid) + '\t' + str(v.length) + '\t' + str(v.cov/(v.length - 63)) + '\n')
# 
# f.close()
# 
# quit()

pbm = PacbioMapping(pbMapping, pbMappingFormat)
print "inputs loaded", len(pbm.readToContig), len(pbm.contigToRead)


# vsSorted = [vid for vid in ig.vs.keys() if vid > 0 and ig.vs[vid].length > 300]
# vsSorted.sort()
# for vid in vsSorted:
#     mc = 0
#     if str(vid) in mv.c2rMapping.readToContig:
#         mc = len(mv.c2rMapping.readToContig[str(vid)])
#     print vid, ig.vs[vid].cov/(ig.vs[vid].length-63), mc
# quit()
#

#
# ig = Graph()
# ig.load('../data/SAureusTW20/SAureusTW20_miseq-contigs.dot',
#         '../data/SAureusTW20/SAureusTW20_miseq-contigs.fa')
# print [(len(v.inne), len(v.oute)) for v in ig.vs.values() if v.vid > 0 and v.length >= 256]
# print len([v for v in ig.vs.values() if v.vid > 0 and v.length >= 256 and str(v.vid) in mv.c2rMapping.readToContig])
# print len([c for c in mv.p2cMapping.contigToRead if ig.vs[int(c)].length >= 256 and len(mv.p2cMapping.contigToRead[c]) > 2])
# print 'Graph loaded', len(ig.vs), len(ig.es)
# pbm = PacbioMapping('../data/SAureusTW20/SAureusTW20_pacbio_miseqcontig_mapping.fasta.m4',
#                     'm4')
# print "inputs loaded", len(pbm.readToContig), len(pbm.contigToRead)
#
#
# ig = Graph()
# ig.load('../data/spalgae/spalgae-contigs.dot',
#         '../data/spalgae/spalgae-contigs.fasta')
# print 'Graph loaded', len(ig.vs), len(ig.es)
# pbm = PacbioMapping('../data/spalgae/sorted_best30.fasta.m4',
#                     'm4algae')
# print "inputs loaded", len(pbm.readToContig), len(pbm.contigToRead)
min_read_threshold = 3
max_read_threshold = 10
nproc = 8
lcg = LongContigGraph(ig, pbm, num_threads=nproc,
                      read_thresholds=(min_read_threshold,
                                       max_read_threshold))
print "LongContigGraph loaded", len(lcg.vs), len(lcg.es)

splcg = lcg.get_simple_paths()
svs = Set([])
svs2 = Set([])
for p in splcg:
    print mv.valid_path(p), p, [lcg.read_count[e.eid] for e in p]
    for v in [p[0].v1, p[-1].v2]:
        if v.Name() in svs:
            svs2.add(v.Name())
    for e in p:
        svs.add(e.v1.Name())
        svs.add(e.v2.Name())
print 'Number of simple paths/ Num valid:', len(splcg),\
len([p for p in splcg if mv.valid_path(p)])
# print splcg
print 'Length 1 paths:', len([s for s in splcg if len(s) == 1])
print 'Num unique contigs covered by simple paths', len(svs)
print 'Total number of contigs in simple paths:',\
sum([len(s) + 1 for s in splcg])
print 'Number of edges/Invalid edges in lcg:', len(lcg.es),\
len([(e, lcg.read_count[e.eid]) for e in lcg.es.values()
     if not mv.valid_path([e])])
print "Invalid edge support",\
[(e, lcg.read_count[e.eid]) for e in lcg.es.values() if not mv.valid_path([e])]
print 'Lengths of path branch vertices:', ([lcg.vs[int(v)].length for v in svs2])
print 'Number of paths with no incoming edge:',\
len([p for p in splcg if len(p[0].v1.inne) == 0])
print 'Number of paths with no outgoing edge:',\
len([p for p in splcg if len(p[-1].v2.oute) == 0])
print 'Paths with exactly 1 incoming edge',\
[(p, p[0].v1.length) for p in splcg if len(p[0].v1.inne) == 1]
print 'Paths with exactly 1 outgoing edge',\
[(p, p[-1].v2.length) for p in splcg if len(p[-1].v2.oute) == 1]

# f = open("../data/spalgae/spalgae_long_contig.dot", 'w')
f = open(outDot, 'w')
# f = open("../data/SAureusTW20/SAureusTW20_long_contig_extended.dot", 'w')
f.write(lcg.__repr__())
f.close()

f = open(outFa, 'w')
# f = open("../data/spalgae/spalgae_long_contig_scaffolds.fa", 'w')
# f = open("../data/SAureusTW20/SAureusTW20_long_contig_extended.fa", 'w')

seen_vertices = {vid: False for vid in ig.vs.keys()}

p_count = max(ig.vs.keys()) + 1
for p in splcg:
    s_start = 0
    s_end = len(p) - 1
#     if len(p[0].v1.inne) > 1 or len(p[0].v1.oute) > 1:
#         s_start = 1
#     if len(p[-1].v2.inne) > 1 or len(p[-1].v2.oute) > 1:
#         s_end = len(p) - 2

    pathseq = ig.path_seq(p[s_start: s_end + 1])
    pathrep = ''
    for e in p[s_start: s_end + 1]:
        pathrep += str(e.v1.vid) + ';' + str(e.ovl) + ';'
    pathrep += str(e.v2.vid)

#     vid = p[0].v1.vid
#     print 'Start', vid, ig.vs[vid].length, vid in lcg.vs,\
#     len(ig.vs[vid].oute),len(ig.vs[vid].inne), len(lcg.vs[vid].oute),\
#     len(lcg.vs[vid].inne)
#     vid = p[-1].v2.vid
#     print 'End', vid, ig.vs[vid].length, vid in lcg.vs,\
#     len(ig.vs[vid].oute),len(ig.vs[vid].inne), len(lcg.vs[vid].oute),\
#     len(lcg.vs[vid].inne)
    f.write('>' + str(p_count) + ' ' + pathrep + '\n')
    f.write(pathseq + '\n') 

    for e in p[s_start: s_end + 1]:
        seen_vertices[e.v1.vid] = True
        seen_vertices[e.v1.conj.vid] = True
        seen_vertices[e.v2.vid] = True
        seen_vertices[e.v2.conj.vid] = True
    p_count += 1

allc = splcg[:]

small_contig_length = 0
for vid in ig.vs.keys():
    if seen_vertices[vid]:
        continue
    if vid < 0:
        continue
    if ig.vs[vid].length < min_vertex_length:
        small_contig_length += ig.vs[vid].length
        continue
#     if vid not in lcg.vs:
#         print vid, ig.vs[vid].length, vid in lcg.vs, len(ig.vs[vid].oute),len(ig.vs[vid].inne), -1, -1, len(mv.c2rMapping.readToContig[str(vid)])
#     else:
#         print vid, ig.vs[vid].length, vid in lcg.vs, len(ig.vs[vid].oute),len(ig.vs[vid].inne), len(lcg.vs[vid].oute),len(lcg.vs[vid].inne), len(mv.c2rMapping.readToContig[str(vid)])
    f.write('>' + str(vid) + '\n')
    f.write(ig.vs[vid].seq + '\n')
    allc.append(ig.vs[vid])
#     print vid
print small_contig_length
f.close()

sev = []
for p in allc:
    if type(p) == list:
        if mv.path_start_dir(p) is None:
            print 'None pathstart', p, len(ig.path_seq(p))
#             if p is None or len(p) <= 3:
#                 continue
#             print p
#             pS, pDir = mv.path_start_dir(p[:3])
#             if pDir == '+':
#                 pE = pS + len(ig.path_seq(p[:3]))
#                 print pS, pE, p[:3]
#             if pDir == '-':
#                 pE = pS - len(ig.path_seq(p[:3]))
#                 print pE, pS, [e.conj for e in p[:3][::-1]]
# 
#             pS, pDir = mv.path_start_dir(p[-3:])
#             if pDir == '+':
#                 pE = pS + len(ig.path_seq(p[-3:]))
#                 print pS, pE, p[-3:]
#             if pDir == '-':
#                 pE = pS - len(ig.path_seq(p[-3:]))
#                 print pE, pS, [e.conj for e in p[-3:][::-1]]

            print [e.conj for e in p[::-1]]
            continue
        pS, pDir = mv.path_start_dir(p)
        if pDir == '+':
            pE = pS + len(ig.path_seq(p))
            sev.append((pS, pE, p))
        else:
            pE = pS - len(ig.path_seq(p))
            sev.append((pE, pS, [e.conj for e in p[::-1]]))
    else:
        print p.vid
        if p.Name() in mv.c2rMapping.readToContig:
            for m in mv.c2rMapping.readToContig[p.Name()]:
                a = PacbioAlignment(m, mv.c2rFileFormat)
                if not a.longAlignment(0.9):
                    continue
                if a.queryStrand == a.seqStrand:
                    sev.append((a.seqStartPos, a.seqEndPos, p))
                else:
                    sev.append((a.seqStartPos, a.seqEndPos, p.conj))

sev.sort(lambda x, y: x[0] - y[0])
print '---------------------------------------------------------------'
print 'Printing sorted sev'
for s in sev:
    if type(s[2]) == list:
        if len(ig.path_seq(s[2])) > min_vertex_length:
            print s
    else:
        print s
print 'sev printing done'
# for si in range(len(sev)-1):
#     olp = sev[si+1][0] - sev[si][1]
#     offset = sev[si][2].length + sev[si+1][0] - sev[si][1]
#     min_offset = 0.9 * offset
#     max_offset = 1.1 * offset
#     print sev[si], sev[si][2].length, sev[si+1][0] - sev[si][1],\
#     ig.find_path(sev[si][2].vid, sev[si+1][2].vid, min_offset, max_offset)

