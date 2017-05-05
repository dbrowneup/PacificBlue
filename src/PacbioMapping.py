#Re-written by Dan Browne on 04/05/16

from datetime import datetime
from math import ceil
import itertools as it
import subprocess as sp
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#from AbstractAlignment import AbstractAlignment



class PacbioMapping:

    def __init__(self, fileName, fileFormat="m4", chunksize=100, num_threads=4):
        print "Entering PacbioMapping module:", str(datetime.now())
        #Read and parse aligns into dictionary of reads (keys) and  alignment objects (values)
        self.flen = self.file_len(fileName)
        self.aligns = iter(open(fileName, "rU"))
        self.format = fileFormat
        self.chunk_limit = chunksize * num_threads
        self.num_chunks = ceil(self.flen / self.chunk_limit)
        self.current_chunk = 0

# Filter reads on the fly in PacbioSubgraph
#        print "Number of unfiltered reads:", len(self.read_iters)
#        print "Number of unfiltered aligns:", len(self.aligns)
#        #Filter out reads with only 1 alignment
##        print "Filtering out reads with only 1 alignment"
#        self.filter_reads()
#        print "Number of filtered reads:", len(self.read_iters)
#        print "Number of filtered aligns:", len(self.aligns)
#        print "Leaving PacbioMapping module:", str(datetime.now())

    def __iter__(self):
        return self

    def next(self):
        if self.current_chunk > self.num_chunks:
            raise StopIteration
        else:
            chunk = dict()
            while len(chunk) <= self.chunk_limit:
                try:
                    A = self.parse_alignment(self.aligns.next().strip('\n').split(' '))
                    try:
                        chunk[A[0]].append(A)
                    except KeyError:
                        chunk[A[0]] = [A]
                except StopIteration:
                    break
            return chunk
        self.current_chunk += 1

#    def filter_reads(self):
#        filtered_reads = set([])
#        for k,v in self.read_iters.items():
#            if len(v) == 1:
#                filtered_reads.add(k)
#                del self.read_iters[k]
#        self.aligns = set([x for x in self.aligns if x.qName not in filtered_reads])

 #   def read_mapping_frequency(self, title, bins=50, range=(0, 50)):
 #       x_axis = [len(x) for x in self.read_iters.itervalues()]
 #       plt.hist(x_axis, bins)
 #       plt.xlabel("aligns per read")
 #       plt.ylabel("Number of reads")
 #       plt.savefig(title+".png")
 #       plt.clf()

    def file_len(self, fname):
        p = sp.Popen(['wc', '-l', fname], stdout=sp.PIPE, stderr=sp.PIPE)
        result, err = p.communicate()
        if p.returncode != 0:
            raise IOError(err)
        return int(result.strip().split()[0])

    def parse_alignment(self, a):
        #Parse alignment into standard tuple:
        #(qName, qStrand, qLength, qStart, qEnd, tName, tStrand, tLength, tStart, tEnd)
        if self.format == 'm4':
            return (str(a[0]), int(a[4]), int(a[7]), int(a[5]), int(a[6]), str(a[1]), int(a[8]), int(a[11]) int(a[9]), int(a[10]))
        elif self.format == 'm5':
            return (str(a[0]), int(a[4]), int(a[1]), int(a[2]), int(a[3]), str(a[5]), int(a[9]), int(a[6]), int(a[7]), int(a[8]))



