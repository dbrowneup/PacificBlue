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
#Contact: kspham@cs.ucsd.edu

from sets import Set

#print args.dotfile
dFile = open('data/spalgae-contigs.dot')

vertexlength = {}
edges = {}

vertexCoverage = {}
overlapLength = {}

for line in dFile:
    if line[0] != '"':
        continue

    if '->' in line:
        v1 = line.strip().split()[0].strip('"')
        v2 = line.strip().split()[2].strip('"')
        if v1 in edges:
            edges[v1].add(v2)
        else:
            edges[v1] = Set([v2])
        d = -63
        if len(line.strip().split()) > 3 and 'd=' in line:
            for i in line.strip().split()[2:]:
                if i.strip('[]')[:2] == 'd=':
                    d = int(i.strip('[]')[2:])
        overlapLength[(v1, v2)] = d

    else:
        v1 = line.strip().split()[0].strip('"')
        ll = line.strip().split()
        for i in ll[0:]:
            if i.strip('[]')[:2] == 'l=':
                vertexlength[v1] = int(i.strip('[]')[2:])
            if i.strip('[]')[:2] == 'C=':
                vertexCoverage[v1] = int(i.strip('[]')[2:])

#ss = 0
#def longNeigbourEdges(v,d):
#    if d == 0:
#        return []
#    founds = []
#    for n in edges[v]:
#        if vertexlength[n] > 1000:
#            founds = founds + [[n]]
#            continue
#        if n not in edges:
#            continue
#        foundxs = longNeigbourEdges(n,d-1)
#        founds = founds + [[n]+i for i in foundxs]
#    return founds

#for v in vertexlength:
#    if vertexlength[v] > 1000 and v in edges:
#        #print v, len(longNeigbourEdges(v,3)), longNeigbourEdges(v,3)
#        ss += len(longNeigbourEdges(v,5))
#
#print sum([len(i) for i in edges.values()]), ss


f = open("data/sorted_best30.fasta.m4")
f.readline()

cread = ""
ccontigs = {}


def concordantContigs(c1, c2):
    for sign1 in "+-":
        if c1[1] + sign1 not in edges:
            continue
        for sign2 in "+-":
            if c2[1] + sign2 in edges[c1[1] + sign1]:
                return True


def constructSubgraph(ccontigs):
    G = {}
    for c1s in ccontigs:
        for c2s in ccontigs:
            if c1s == c2s:
                continue
            c1 = c1s.strip().split()
            c2 = c2s.strip().split()
            if concordantContigs(c1, c2):
                if c1[1] in G:
                    G[c1[1]].add(c2[1])
                else:
                    G[c1[1]] = Set([c2[1]])
    return G


def isLineGraph(G):
    if len(G) == 0:
        return False
    c0 = G.keys()[0]
    q = [c0]
    qi = 0
    while qi < len(q):
        #print q, qi, G[q[qi]]
        #if q[qi] not in G: continue
        for c in G[q[qi]]:
            if c not in q:
                q.append(c)
        qi += 1
    if len(q) == len(G.values()) and \
       len([c for c in G if len(G[c]) == 1]) == 2 and \
       len([c for c in G if len(G[c]) > 2]) == 0:
        return True
    else:
        return False


subgraphList = {}
isLineList = {}
contigstoPacbio = {}  # [contig]  Set("lines")

for line in f:
    ll = line.strip().split()
    if ll[1] in contigstoPacbio:
        contigstoPacbio[ll[1]].add(line)
    else:
        contigstoPacbio[ll[1]] = Set([line])
    if ll[0] != cread:
        subgraphList[cread] = constructSubgraph(ccontigs)
        isLineList[cread] = 0
        if isLineGraph(subgraphList[cread]):
            isLineList[cread] = len(subgraphList[cread])
        cread = ll[0]
        ccontigs = Set([])
    ccontigs.add(line)

subgraphList[cread] = constructSubgraph(ccontigs)


def VLengthCompare(v1, v2):
    return vertexlength[v1] - vertexlength[v2]


def getCandidateSeeds(vertexlength):
    V = [v for v in vertexlength.keys() if vertexlength[v] > 3000]
    V.sort(VLengthCompare)
    return V


candidates_seeds = getCandidateSeeds(vertexlength)


def Extend(seed, avoid):
    if seed not in edges.keys():
        return [seed]
    else:
        print "Im extending", seed
    if seed[-1] == "+":
        maplist = []
        for mapping in contigstoPacbio[seed[:-1]]:
            rr = mapping.strip().split()
            if rr[8] == '1' and int(rr[7]) - int(rr[6]) > 100:
                maplist.append(mapping)
            elif rr[8] == '0' and int(rr[5]) > 100:
                maplist.append(mapping)
        neiborlist = Set([])
        for r in maplist:
            for c in subgraphList[r.strip().split()[0]]:
                if c == seed:
                    continue
                if c + '+' in edges[seed]:
                    neiborlist.add(c + '+')
                if c + '-' in edges[seed]:
                    neiborlist.add(c + '-')
        if len(neiborlist.difference(avoid)) == 0 or \
           len(neiborlist.difference(avoid)) > 1:
            print "I can not find any extension", seed
            return [seed]
        else:
            print "I can find one extension for", seed, (list(neiborlist))[0]
            return [seed] + \
                Extend((list(neiborlist.difference(avoid)))[0],
                       Set([(list(neiborlist.difference(avoid)))[0]]
                       + list(avoid)))
    return [seed]

covered_seeds = Set([])
for seed in candidates_seeds:
    if seed in covered_seeds:
        continue
    contigForward = Extend(seed, Set([seed]))
    covered_seeds.union(Set(contigForward))
