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
#Contact: kspham@cs.ucsd.edu

import logging


class Abstract_Vertex(object):
    def __init__(self, vid):
        self.inne = []  # incoming edges
        self.oute = []  # outgoing edges
        self.vid = vid  # vertexid

    def outv(self):
        return [e.v2 for e in self.oute]

    def innv(self):
        return [e.v1 for e in self.inne]

    def __hash__(self):
        return self.vid


class Abstract_Edge(object):
    def __init__(self, v1, v2, eid):
        self.v1, self.v2 = v1, v2
        v1.oute.append(self)
        v2.inne.append(self)
        self.eid = eid

    def __hash__(self):
        return self.eid

    def length(self):
        pass


class Abstract_Graph(object):

    def __init__(self):
        self.es = {}  # key -->edges
        self.vs = {}  # key -->vertices
        self.logger = logging.getLogger('Algae')
