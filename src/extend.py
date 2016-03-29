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
import operator
from PacbioMapping import PacbioMapping
from PacbioAlignment import PacbioAlignment
from PacbioSubgraph import PacbioSubgraph
from loadgraph import loadgraph
from illumina_graph import Graph
from validation.mapping_validation import MappingValidation

g = Graph()
g.load(experimental.illumina_graph_file, experimental.illumina_contig_file)
mv = MappingValidation(experimental.contig_refmap_file, experimental.contig_pacbiomap_file, experimental.pacbio_refmap_file)
g.load_dna(experimental.illumina_contig_file)
m = PacbioMapping(experimental.contig_pacbiomap_file)
ma = PacbioMapping(experimental.contig_refmap_file, 'blastn')

idlst = {}
visited = []
overlap = {}
len_count = 0

for vid, v in g.vs.iteritems(): 
    if v.length > 5000:
        len_count += 1
    idlst[vid] = v.length

sorted_idlst = sorted(idlst.iteritems(), 
                      key=operator.itemgetter(1), 
                      reverse=True)

#print sorted_idlst
#for x in  sorted_idlst:
#    print x
#    for y in g.Subgraph(x[0], 1):   
#        for z in y:
#            if z.ovl < -500 and z.v1 not in visited:
#                visited.append(z.v1)
#                overlap[z.v1.vid] = z.ovl
#sorted_overlap = sorted(overlap.iteritems(), 
#                        key=operator.itemgetter(1), 
#                        reverse=False)
#print sorted_overlap

def ck_intersect(paID, contig_id, key = None):
#    print paID
    pac_subg = PacbioSubgraph(paID, g, m, augment= False)    
    if key is None:
#        if len(pac_subg.es) < 3 and len(pac_subg.es) > 0:
#    #        print paID
##            print pac_subg.vs.values()
#            for outpaths in pac_subg.es:
##                print len(pac_subg.es)                                           
#       #     if len(pac_subg.es) < 3 and len(pac_subg.es) > 0 :
#                if abs(pac_subg.es[outpaths].v2.vid) == abs(contig_id):
#                    return -1 * pac_subg.es[outpaths].v1.vid
#                elif abs(pac_subg.es[outpaths].v1.vid) == abs(contig_id):
#                    return pac_subg.es[outpaths].v2.vid
#    else:       
#
#        if len(pac_subg.es) > 3 and len(pac_subg.es) < 5:
#            a = []
#            for path in pac_subg.es:
#                next_contig = 0
#                if abs(pac_subg.es[path].v1.vid) != abs(contig_id):
#                    a.append(pac_subg.es[path].v1.vid)
#    #        print a
#            for x in a:
#                if (-1 * x) not in a:
#    #                print  x
#                    return x
        if len(pac_subg.es) > 0:
            a = []
            b = []
#            print pac_subg.es
            for paths in pac_subg.es:
                if pac_subg.es[paths].v1.vid == contig_id:
                    a.append(pac_subg.es[paths].v2.vid)
                    b.append(pac_subg.es[paths].v2.vid)
#                    print pac_subg.es[paths]
            for x in a:
                for paths in pac_subg.es:
                    if pac_subg.es[paths].v1.vid == x:
                        b.append(pac_subg.es[paths].v2.vid)
#                        print pac_subg.es[paths]
            return b    



def pacbio_extend(contig_id, key = None):
    graphs = []
    vertices = {}
#    print m.contigToRead.keys()
    if str(abs(contig_id)) in m.contigToRead.keys():
        for i in xrange(len(m.contigToRead[str(abs(contig_id))])):
            pa = PacbioAlignment(m.contigToRead[str(abs(contig_id))][i],'m4')
            paID = pa.queryID
            ck_pac_reads = ck_intersect(paID, contig_id, key)
            if ck_pac_reads != None and len(ck_pac_reads) > 0:
                for y in ck_pac_reads:
                    graphs.append(y)
        
        for x in graphs: 
            if abs(x) != abs(contig_id):
                if x not in vertices:
                    vertices[x] = 0
                vertices[x] += 1
#        for gg in graphs:
#            for i, vid in enumerate(gg.vs.keys()):
#                if abs(vid) != abs(contig_id) and check_overlap(contig_id, vid):
#                    if vid not in vertices: 
#                        vertices[(vid)] = 0
#                    vertices[(vid)] += 1
##                    print pacid[i], vid
#        sorted_vertices = sorted(vertices.iteritems(), 
#3                                 key=operator.itemgetter(1), 
#3                                 reverse=True)
#        return sorted_vertices
    return vertices
    

#print pacbio_extend(664)
#print pacbio_extend(-657,-668)
def check_overlap(contig_id, out_contig_id):
    for outpaths in g.Subgraph(contig_id, 1):
        for path in  outpaths:
            if path.v2.vid == out_contig_id and path.ovl < -200:
                return True 
                
    return False            


#Extend the outgoing paths from a given contig id and the 
#radius that wish to travel. Output as Dot format
def make_graph_from_contig(contig_id, radius):
    outgoing_path = []
    if len(g.Subgraph(contig_id, radius)) > 0: 
        with open(experimental.path_extend_graph_file, 'w') as graph_file:    
    #        if len(g.Subgraph(contig_id, radius)[0]) > 0:
            graph_file.write('digraph paths {')
#            print g.Subgraph(contig_id,radius)
            graph_file.write('graph [k=%s]' % \
                            (len(g.Subgraph(contig_id,radius)[0])))
       
            vertices = [] 
            for i, outpaths in enumerate(g.Subgraph(contig_id,radius)):
                for path in outpaths:
                    if i == 0: 
                        outgoing_path.append(path.v2.vid)
    #                if path.v1.vid == -657 and path.v2.vid == -668:
    #                    print path.v1.vid, path.v2.vid, path_color(ma, [path])
                    graph_file.write('"%s%s" -> "%s%s" [color=%s] [label=%s]' % \
                                    (path.v1.Name(), path.v1.Strand(),
                                     path.v2.Name(), path.v2.Strand(),
                                     path_color(ma, [path]),
                                     path.ovl))
                                     
                    if path.v1 not in vertices:
                        vertices.append(path.v1)
                    if path.v2 not in vertices:
                        vertices.append(path.v2)
    
            for vertex in vertices:
                graph_file.write('"%s%s" [label="%s%s l=%s "]' % \
                                (vertex.Name(), vertex.Strand(),
                                 vertex.Strand(), vertex.Name(),
                                 vertex.length))
#                                 ma.positions_on_seq(vertex.Name(),
#                                                     vertex.Strand())))
            graph_file.write('}')
        graph_file.close()
    return outgoing_path 


def path_color(ma, path):
    if mv.valid_path(path):
        return "red"
    return "black"

def get_out_id(contig_id):
#    print contig_id
    outid = []
     
    if contig_id in g.vs.keys():
        
        for x in g.vs[contig_id].oute:
#        print x
            outid.append(x.v2.vid)
        return outid
    else:
        return None

def pick_path(contig_id, pre_contig = None):
    
    ig = get_out_id(contig_id)
    if ig == None:
        return None

    else:
        pg = []
        pg_dict = pacbio_extend(contig_id)
        pg = pg_dict.keys()
    
        pg = set(pg) 
        ig_pg_intersect = list(set(ig) & pg)
#        print contig_id
#        print ig
#        print pg_dict
#        print ig_pg_intersect
        if pre_contig is not None:
            pre_pg_dict = pacbio_extend(pre_contig)
            for k in pre_pg_dict.keys():
                if k in pg:
                    pg_dict[k] += pre_pg_dict[k]
                
        if len(ig_pg_intersect) > 0:
            maxid = 0
            maxcount = -1

            sec_id = 0
            sec_count = -1
           
            for x in ig_pg_intersect:  
                if x in pg_dict.keys():
                    if len(ig_pg_intersect) > 1:
                        if pg_dict[x] > sec_count:            
                            sec_count = pg_dict[x]
                            sec_id = x
                            if sec_count > maxcount:
                                sec_count = maxcount
                                sec_id = maxid
                                maxcount = pg_dict[x] 
                                maxid = x
                    else:
                        if pg_dict[x] > maxcount:
                            maxcount = pg_dict[x]
                            maxid = x
            
#            print maxcount, sec_count, "max and sec"                    
            if maxcount - sec_count  <= 1 :
#                print maxid
                return None
#                print maxid, " merge",'\n'
            elif maxcount > 3:
                return maxid

            elif (pg_dict.keys()[pg_dict.values().index(max(pg_dict.values()))] in 
                  get_out_id(maxid)):
                if max(pg_dict.values()) > 2:
                    return maxid

            elif pre_contig is not None: 
                return double_ck(maxid, pre_contig) 

            else:
                return None
        else:
#            print ig
#            if (len(ig) == 1 and 
#                g.vs[ig[0]].length < 700 and 
#                g.vs[ig[0]].length > 550 and 
#                g.vs[contig_id].length > 100 and
#                g.vs[ig[0]].length - g.vs[contig_id].oute[0].ovl > 20 ):
#                return "ig[0]"
                
            if pre_contig is not None:                
                return double_ck(contig_id, pre_contig)
            return None
        

#print pick_path(-432)
def double_ck(contig, pre_contig):
    
    intersect = []
#   if pre_contig is not None:
    contig_dict = pacbio_extend(contig)
    pre_contig_dict = pacbio_extend(pre_contig)
#    print pre_contig , contig," pre" , "contig"
#    print contig_dict," contig"
#    print (pre_contig_dict),"pre" + "\n"

    for contigs in contig_dict.keys():
        if contigs in pre_contig_dict.keys():
#            print contigs, " contig"
            intersect.append(contigs)


#    print intersect ,"intersect"
    if len(intersect) > 0 :
        maxcount = 2 
        maxinter_id = 0
#        print intersect," intersect"
        the_count = []
        
        for x in intersect:
#            print pre_contig_dict[x], contig_dict[x]," the count"
            reads_count = int(pre_contig_dict[x]) - int(contig_dict[x])
            the_count.append(reads_count)
            if reads_count > maxcount:
                maxcount = reads_count
                maxinter_id = x
#        print maxinter_id, "inter id"
        if maxinter_id != 0:
            return maxinter_id

#    if len(intersect) == 0 or pre_contig is not None: 
#        print "work", pre_contig, contig
#        ig = []
#
#        for x in g.Subgraph(contig, 5):
#            for y in x:
#                ig.append(y.v2.vid)
#
#        print ig 
#        extra_reads = pacbio_extend(pre_contig, key="double check")
#        for k in extra_reads.keys():
#            if k in ig:
#                return k 
    return None
#print double_ck(-743, 741)
#print pick_path(-724)
def search_out_path(contig):
    temp = []
    def dfs(contig): 
        if len(temp) > 1 :
            next_contig = pick_path(contig, temp[len(temp)-1])
        else:
            next_contig = pick_path(contig)

        if (next_contig == None):
#        or g.vs[next_contig].length < 300):
       # or 
        #    next_contig in visited):
#            if len(temp) > 1:
#                if temp[len(temp)-1] in visited:
#                    visited.remove(temp[len(temp)-1])
#                elif (-1 * temp[len(temp)-1]) in visited:
#                    visited.remove(-1 * temp[len(temp)-1])
#
            return 
        else:
            if (next_contig != None and 
                next_contig != 0 and 
                next_contig not in visited):
                visited.append(next_contig)
                visited.append(-1 * next_contig)
                temp.append(next_contig)
                dfs(next_contig)
    temp.append(contig)
#    visited.append(contig)
    dfs(contig)
#    print contig, visited
#    if len(temp) == 1:
#        visited.remove(temp[0])
#        visited.remove(-1*temp[0])
    return temp

#print search_out_path(720)
def select_from_pac(contig):
    
    out = pacbio_extend(contig)
#    print out 
    if len(out) == 0:
        return None
    elif len(out) == 1 and out.values()[0] > 5:
        return out.keys()[0]
    else: 
        return None
    
#print a 
    if len(out) > 1:
        a = sorted(out.values(), reverse=True)
#        print a 
        if a[0] - a[1] > 3:
            return out.keys()[out.values().index(a[0])]
        else:
            return None

def search_out_pac(contig):
    temp = []
    
    next_c = select_from_pac(contig)
#    if next_c != None:
#        temp.append(contig)
    while next_c != None:
        temp.append(next_c)
        next_c = select_from_pac(next_c)
    return temp

#print pacbio_extend(350)
#print search_out_pac(657)
#print select_from_pac(657)
#print search_out_path(227)," the path"
#print pick_path(-668, -657)," pick"
def run(contigs_list):
#    print len(contigs_list), "lenght of contigs_list"
    contig_ext = []
#    boolst = [False] * len(contigs_list) 
#    visited = []
    for i in xrange(len_count-1):
#    for i in xrange(len(contigs_list)):
        if contigs_list[i] not in visited or abs(contigs_list[i] not in visited):
#            visited.append(contigs_list[i])
            
            full_path = search_out_path(contigs_list[i])
#            print visited 
            if len(full_path) > 0 and full_path not in contig_ext:
                contig_ext.append(full_path)


    for k in xrange(len_count, len(contigs_list)-1):
        if contigs_list[k] not in visited:
            contig_ext.append([contigs_list[k]])
#        print len(visited), "length of visited"
#        print len(contig_ext), " length of contig_ext"
    return contig_ext


def count_read_map(contig_id):
    temp = {}
    reads = []
    if str(abs(contig_id)) in m.contigToRead.keys():
        for i in xrange(len(m.contigToRead[str(abs(contig_id))])):
            pa = PacbioAlignment(m.contigToRead[str(abs(contig_id))][i],'m4')
            paID = pa.queryID
#            print paID, "IDDD"
             

            for j in xrange(len(m.readToContig[paID])):
                pa1 = PacbioAlignment(m.readToContig[paID][j], 'm4')
                
                a = pa1.seqStrand + pa1.seqID
                if int(a) == contig_id and pa1 not in reads:
                    reads.append(pa1)
                
           
#        print len(reads), "read"
        for x in reads:
            for k in xrange(len(m.readToContig[x.queryID])): 
                pa2 = PacbioAlignment(m.readToContig[x.queryID][k], 'm4')
                b = pa2.seqStrand + pa2.seqID
                if b not in temp:
                    temp[b] = 0
                temp[b] += 1
   
#        c = sorted(temp.items(), key=operator.itemgetter(1), reverse = True)
#        print c
        return temp
        

#        print b

                
#print get_out_id(723) 
#count_read_map(-680)
def get_offset(c1, c2):
    a1 = count_reads_map(c1)
    a2 = count_reads_map(c2)
    c = list(set(a1) & set(a2))
    a1_alig = get_alig(c1, c)
    a2_alig = get_alig(c2, c)
    offset = []
    for j in c:
        
        for k in a1_alig[j]:
            for l in a2_alig[j]:
                d =  k.relative_offset(l, g.vs[c1], 
                                        g.vs[c2])
                if d != False:
                    offset.append(d)
    if len(offset) > 2 :
        offset_ave = reduce(lambda x, y: x + y, offset) / len(offset)
    
        t_offset = offset_ave - g.vs[c1].length
        return t_offset


def extend_gap(contigs):
    temp = []
    contig = contigs[len(contigs)-1]
#    print contig
    a = sorted(count_read_map(contig).iteritems(), key = operator.itemgetter(1), reverse = True)
#    print a 
    for x in a:
        
#        print type(x[0]), type(contig), type(x[1])
        if abs(int(x[0])) != abs(int(contig)) and int(x[1]) > 0: 
            temp.append(x)
        
    #print temp, "temp"
#    abc = []
#    for x in temp:
#       abc.append(abs(int(x[0])) in visited)
#    print abc, " abc"
#    print visited

#    for x in temp:
#        if int(x[0]) not in visited:
#            if contig < 0:
#                return -1 * x[0]
#            else:
#                return x[0]
    for y in temp:
        if int(y[0]) in contigs:
            temp.remove(y)
#    if contig == -698:
#    print temp

#    for k in temp:
#        inne = []
#        oute = [] 
#           
#        for z in g.vs[int(k[0])].inne:
#            inne.append(z.v1.vid)
#        for s in g.vs[int(k[0])].oute:
#            oute.append(s.v2.vid)
##       
#        print (k[0], g.vs[int(k[0])].length, 
#                int(contig) in (inne), 
#                int(contig) in (oute))
###        if int(contig) in (inne):
##            return None
    if len(temp) > 3: 
        if (temp[0][1] - temp[1][1] > 1):
       # and temp[0][1] < 14):
       # and int(temp[0][0]) not in visited:
            offset = get_offset(int(contig), int(temp[0][0]))
            if ((offset > 30) or
               (offset < -50 and 
                offset > -5000)):

                return temp[0][0]

#            if int(contig) < 0:
#                return -1 *int(temp[0][0])
#            else:
#                return temp[0][0]
#            if int(contig) < 0:
#                return -1 *int(temp[0][0])
#        else:
#            minid = -1
#            mincount = 999
##            print temp
#            for y in temp:
#                if len(g.vs[int(y[0])].inne) < len(g.vs[int(y[0])].oute):
#                    count = len(g.vs[int(y[0])].inne) + len(g.vs[int(y[0])].oute)
#                    if count < mincount:
#
#                        minid = y[0]
#                        mincount = count
##                        print y[0], len(g.vs[int(y[0])].inne), len(g.vs[int(y[0])].oute), g.vs[int(y[0])].inne, g.vs[int(y[0])].oute 
#            if mincount < 4:
#                return minid 
#            return None
    return None
    
#    print temp
   # break
    

#    print a
#print g.vs[734].length 
#print extend_gap(710)

def fill_gaps(merged_contigs):
    head_list = []
    tail_list = []
    
    for k in merged_contigs:
        if len(k) > 1:
            head_list.append(k[0])
            tail_list.append(k[len(k)-1])
    temp = [] 
    for i, x in enumerate(merged_contigs):  
        if len(x) > 1:
            a =  search_out_pac(-1 * x[0])
            if len(a) > 1:
                for j in a:
                    x.insert(0, -1*j)

            b = search_out_pac(x[len(x)-1])
            if len(b) > 1:
                for k in b:
                    x.append(k)

        elif len(x) == 1:    
            c = search_out_pac(x[0])
            if len(c) > 1:
                for p in c:
                    x.append(p)
#        print type(x[0]), x
        if len(x) > 1:
#            next_head = extend_gap(x[0])
            next_tail = extend_gap(x)
#            while next_head != None:
#                if next_head != "":            
#                    
#                    x.insert(0, (-1 * int(next_head)))
#                    print x
#                    
#                    next_head = extend_gap(x[0])
#    
#                else:
#                    break
            added = []
            while next_tail != None:
                if (next_tail != "" and next_tail not in added):
                    added.append(next_tail)

                    x.append(int(next_tail))
#                    print x
                    next_tail = extend_gap(x)
                else:
                    break
        
        temp.append(x)
    return temp
#print extend_gap(-750)           
def count_reads_map(contig_id):
    temp = []
    if str(abs(contig_id)) in m.contigToRead.keys():
        for i in xrange(len(m.contigToRead[str(abs(contig_id))])):
            pa = PacbioAlignment(m.contigToRead[str(abs(contig_id))][i],'m4')
#             print pa
            paID = pa.queryID
             
            if paID not in temp:
                temp.append(paID)
    return temp


def get_alig(contig_id, reads):
    temp = {}
    alig = []
    if str(abs(contig_id)) in m.contigToRead.keys():
        for i in xrange(len(m.contigToRead[str(abs(contig_id))])):
            pa = PacbioAlignment(m.contigToRead[str(abs(contig_id))][i],'m4')
            if pa.queryID in reads:
                if pa.queryID not in temp.keys():                    
                    temp[pa.queryID] = []
                temp[pa.queryID].append(pa)
           
#            if paID not in temp:
#               temp.append(paID)
    return temp
  
#print get_offset(521, 734)
#print pacbio_extend(749)
#make_graph_from_contig(-698, 5)
def get_edges(contig_ids):
    edges = []
    
    if len(contig_ids) == 1:
        return contig_ids[0]
    for i in xrange(len(contig_ids)-1):
        added = False
        for contig in g.vs[contig_ids[i]].oute:
            if (contig_ids[i] == contig.v1.vid and
                contig_ids[i + 1] == contig.v2.vid):
                edges.append(contig)
                added = True
        if added == False:
#            print (" please adddddddddddddddddddddddddddddddd")
#            a = []
            a1 = count_reads_map(contig_ids[i])
            a2 = count_reads_map(contig_ids[i+1])
            c = []
#            for x in a1:
#                if x in a2:
#                    c.append(x)
#            print c, "ccccccc"
#            print a1
#            print a2
            c = list(set(a1) & set(a2))
#            print c , "cccc"
#            print contig_ids[i], contig_ids[i+1]
            a1_alig = get_alig(contig_ids[i], c)
            a2_alig = get_alig(contig_ids[i+1], c)
#            print a1_alig[c[2]]
#            print a2_alig[c[2]] 
            offset = []
            offsetc = []
            for j in c:
                
                for k in a1_alig[j]:
#                    print a2.alig[j]
                    for l in a2_alig[j]:
#                        print l, "LLLLLL"
                        d =  k.relative_offset(l, g.vs[contig_ids[i]], 
                                                g.vs[contig_ids[i+1]])
                        conj = l.relative_offset(k, g.vs[contig_ids[i+1]].conj,
                                                 g.vs[contig_ids[i]].conj)


#                        print d
                        if d != False:
                            offset.append(d)
                        if conj != False:
                            offsetc.append(conj)
            if len(offset) > 2 and len(offsetc) > 2:
                offset_ave = reduce(lambda x, y: x + y, offset) / len(offset)
                offsetc_ave = reduce(lambda x, y: x + y, offsetc) / len(offsetc)
    
                t_offset = offset_ave - g.vs[contig_ids[i]].length
                t_offsetc = offsetc_ave - g.vs[contig_ids[i+1]].length
    
                g.add_edge(-1*contig_ids[i+1],-1*  contig_ids[i], t_offsetc)
                g.add_edge(contig_ids[i], contig_ids[i+1], t_offset)
                edges.append(g.add_edge(contig_ids[i], contig_ids[i+1], t_offset))
    return edges

#tlist = [664, -728] 
#print get_edges(tlist)

def make_fasta(sorted_idlst):
    
    idonly = []
    for y in sorted_idlst:
        idonly.append(y[0])
    
    pre_merged_contig = run(idonly)
    
#    for x in pre_merged_contig:
#        print x
#    print pre_merged_contig
    def merged(pre_merged):
        merged_contig = []
#        for i in xrange(len(pre_merged) - 2):
#            temp = []
#            if (abs(pre_merged[i][0]) == abs(pre_merged[i+1][0])):
#        #        print pre_merged_contig[i]
#                a = pre_merged[i]
#                if len(a) > 1:
#                    for j in xrange(len(a)):
#                        temp.append(a[j] * -1)
#                else:
#                    for j in a:
#                        temp.append(j)
#                    
#                c = temp[::-1] + pre_merged[i+1][1::]
#                merged_contig.append(c)
#            else:
#                added.append(i) 
                
    #    print len(merged_contig)
#        added = []
        posit = []
        negat = []
        for k in pre_merged:
            if k[0] > 0:
                posit.append(k)
            else:
                negat.append(k)
        for p in posit:
            for n in negat:
                if p[0] == abs(n[0]):                    
                    temp = []
                    a = p
                    for j in xrange(len(a)):
                        temp.append(a[j] * -1)
                    c = temp[::-1] + n[1::]
                    merged_contig.append(c)
                    posit.remove(p)
                    negat.remove(n)
        for x in posit:
            merged_contig.append(x)
        for y in negat:
            merged_contig.append(y)
#        print len(merged_contig), len(posit), len(negat)
        return merged_contig
    
#    for x in pre_merged_contig:
#        print x
    merged_contig = merged(merged(pre_merged_contig))
#    print merged_contig
    filled_merged_contig = fill_gaps(merged_contig)
#    print filled_merged_contig
#    merged_edges = []
#    
    with open("jointed_contig.fasta", "w") as output:
        for x in filled_merged_contig:
            edges = get_edges(x)
            if len(x) > 1 :
#                print x
#                print mv.valid_path(edges), edges
                output.write("> %s\n" % x)
                output.write("%s\n" % g.path_seq(edges))   
                
#                merged_edges.append(edges)
#                print mv.valid_path(edges), edges
#                print edges
#                print x 
            else:
#                print x , g.vs[x[0]].length
                if x[0] not in visited or abs(x[0]) not in visited:
#                and g.vs[x[0]].length > 2000:
                    
#                    print x , g.vs[x[0]].length
                    output.write("> %s\n" % x)
                    output.write("%s\n" % g.vs[x[0]].seq)
##print sorted_idlst 
make_fasta(sorted_idlst)


#print extend_gap([-657,-668, 664])





#make_graph_from_contig(749, 4)
#print g.vs[683].length, "683"
#print g.vs[746].length, "746" 
#print g.vs[670].length, "670"


#print len(merged_contig)
#test =[[-712, -688, -673]]
#with open('joined_contig.fasta', 'w') as output:
#    for x in merged_contig:       
##        if len(x) > 1:
##        print '> %s\n' % (x)
##        print '%s\n' % str(g.get_seq(x))
#        output.write('> %s\n' % (x))  
#        output.write('%s\n' % str(g.get_seq(x)))
#        print x
'''
premerged = merged(idonly)
with open("before_joined.fasta", 'w') as output1:
    for x in premerged:       
#       if len(x) > 1:
       output1.write('> %s\n' % (x))  
       output1.write('%s\n' % str(g.get_seq(x)))
       
'''
   

#print g.vs[-1].seq
#print g.vs[1].seq
#print g.vs[1].seq[::-1]
#print g.vs[1].seq[::-1][3::]









#
#vert_to_red = [521, 19,179,619,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,670,672,673,674,675,676,677,678,679,680,681,682,683,684,686,687,688,689,690,691,692,693,694,695,696,697,698,699,701,702,703,705,708,709,710,711,712,713,714,716,717,719,721,723,724,725,726,728,731,732,733,734,735,736,737,738,740,741,743,744,745,746,747,748,749,750,751,752,-752,-751,-750,-749,-748,-747,-746,-745,-744,-743,-741,-740,-738,-737,-736,-735,-734,-733,-732,-731,-728,-726,-725,-724,-723,-721,-719,-717,-716,-714,-713,-712,-711,-710,-709,-708,-705,-703,-702,-701,-699,-698,-697,-696,-695,-694,-693,-692,-691,-690,-689,-688,-687,-686,-684,-683,-682,-681,-680,-679,-678,-677,-676,-675,-674,-673,-672,-670,-668,-667,-666,-665,-664,-663,-662,-661,-660,-659,-658,-657,-656,-655,-654,-653,-619,-179,-19,-521]
#
#
#
##contig_ids = [-619, -733, 661, -663, -673, -674, 687, -675, 703, -677, 694, 691, 746, 747, -684, 726, 702, -708, -709, -710, 717, -723, 746, -716, -719, -724, 734, 655, 748, 751, -752, 692, 726, 667]
#
#
#def make_graph(contig_ids, radius):
#    outgoing_path = []
#    for contig_id in contig_ids:
#        if len(g.Subgraph(contig_id, radius)) > 0: 
#            with open("dead_end_graphs/C" + str(contig_id) + ".dot", 'w') as graph_file:    
#        #        if len(g.Subgraph(contig_id, radius)[0]) > 0:
#                graph_file.write('digraph paths {')
#    #            print g.Subgraph(contig_id,radius)
#                graph_file.write('graph [k=%s]' % \
#                                (len(g.Subgraph(contig_id,radius)[0])))
#           
#                vertices = [] 
#                for i, outpaths in enumerate(g.Subgraph(contig_id,radius)):
#                    for path in outpaths:
#                        if i == 0: 
#                            outgoing_path.append(path.v2.vid)
#        #                if path.v1.vid == -657 and path.v2.vid == -668:
#        #                    print path.v1.vid, path.v2.vid, path_color(ma, [path])
#                        graph_file.write('"%s%s" -> "%s%s" [color=%s] [label=%s]' % \
#                                        (path.v1.Name(), path.v1.Strand(),
#                                         path.v2.Name(), path.v2.Strand(),
#                                         path_color(ma, [path]),
#                                         path.ovl))
#                                         
#                        if path.v1 not in vertices:
#                            vertices.append(path.v1)
#                        if path.v2 not in vertices:
#                            vertices.append(path.v2)
#        
#                for vertex in vertices:
#                    if vertex.vid in vert_to_red:
#                        color = "red"
#                    else:
#                        color = "black"
#                    graph_file.write('"%s%s" [label="%s%s l=%s"]' % \
#                                    (vertex.Name(), vertex.Strand(),
#                                     vertex.Strand(), vertex.Name(),
#                                     vertex.length ))
#
#                graph_file.write('}')
#            graph_file.close()
#        print str(contig_id) + "done"
#    return outgoing_path 
#
#
#make_graph(contig_ids, 4)
#
#

