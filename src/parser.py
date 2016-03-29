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

#Author: Son Pham
#Contact: kspham@eng.ucsd.edu

from utils import rc
from path.path import path

#parse graph vertices


def graph_vertices(filename):
    filepath = path(filename)
    graph = filepath.open()
    ignorelines = 3
    for _ in range(ignorelines):
        graph.readline()

    while True:
        line = graph.readline()
        graph.readline()
        if '->' in line:
            break
        vinf = line.strip().split()
        v = vinf[0].strip('"')
        vid = int(v[:-1])
        if vid == 0:
            vid = 721199353
        if v[-1] == '-':
            vid = -1 * vid
        length = int(vinf[1][3:])
        cov = int(vinf[2][2:-1])
        vid_conj = -1 * vid
        yield vid, vid_conj, length, cov

    graph.close()

#parse graph edges


def graph_edges(filename):
    filepath = path(filename)
    graph = filepath.open()
    for _ in range(3):
        line = graph.readline()
    defaultoverlap = int(line.strip().split()[1][3:-1])

    while '->' not in line:
        line = graph.readline()
    while True:
        edgeinfo = line.strip().split()
        v1id = int(edgeinfo[0][1:-2])
        if v1id == 0:
            v1id = 721199353
        if edgeinfo[0][-2] == '-':
            v1id = -1 * v1id
        v2id = int(edgeinfo[2][1:-2])
        if v2id == 0:
            v2id = 721199353
        if edgeinfo[2][-2] == '-':
            v2id = -1 * v2id

        newoverlap = defaultoverlap
        if len(edgeinfo) == 4:
            newoverlap = int(edgeinfo[3][3:-1])
        yield v1id, v2id, newoverlap
        line = graph.readline()
        if '}' in line:
            break
    graph.close()


def contigs_sequence(filename):
    filepath = path(filename)
    contigs = filepath.open()
    while True:
        idinf = contigs.readline().strip()
        if len(idinf) < 1:
            break
        seq = contigs.readline().strip()
        contigid = int(idinf.strip().split()[0][1:])
        if contigid == 0:
            contigid = 721199353
        yield contigid, seq
        yield -1 * contigid, rc(seq)
