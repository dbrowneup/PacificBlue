#Re-written by Dan Browne on 04/05/16

from AbstractAlignment import AbstractAlignment
from datetime import datetime
import matplotlib
import matplotlib.pyplot as plt

class PacbioMapping:

    def __init__(self, fileName, fileFormat="m4"):
        print "Beginning to load PacBio mapping:", str(datetime.now())
        self.readToContig = {}
        self.alignments = open(fileName, "rU").read().split('\n') # reads alignments
        self.alignments = self.alignments[1:-1] if "score" in self.alignments[0] else self.alignments[0:-1] # discards empty last row and first row if header is detected
        self.alignments = [AbstractAlignment(x.split(' '), fileFormat) for x in self.alignments] # parses alignments into AbstractAlignment objects
        self.readToContig = {x.qName: [] for x in self.alignments if x.qName not in self.readToContig} # parses qNames into dictionary
        #parse alignments into dictionary
        for x in self.alignments:
            self.readToContig[x.qName].append(x)
        print "Number of mapping reads:", len(self.readToContig)
        print "Number of alignments:", len(self.alignments)
        #sort readToContig mappings in decreasing order of alignment length
        for k in self.readToContig:
            self.readToContig[k] = sorted(self.readToContig[k], key=lambda x: (x.qEnd - x.qStart), reverse=True)
        print "Finished loading PacBio mapping:", str(datetime.now())
    def filter_reads(self, length_fraction=0.7):
        #remove all reads with only 1 alignment
        #remove reads where 1st alignment covers >= length_fraction (default 70%) of read length
        print "Beginning read filtration:", str(datetime.now())
        filtered_reads = set([])
        for k,v in self.readToContig.items():
            first_map_length = v[0].qEnd - v[0].qStart
            if len(v) == 1:
                filtered_reads.add(k)
                del self.readToContig[k]
            elif first_map_length >= v[0].qLength * float(length_fraction):
                filtered_reads.add(k)
                del self.readToContig[k]
        self.alignments = set([x for x in self.alignments if x.qName not in filtered_reads])
        print "Number of mapping reads filtered out:", len(filtered_reads)
        print "Number of mapping reads remaining:", len(self.readToContig)
        print "Finished read filtration:", str(datetime.now())
    def read_mapping_frequency(self, title, bins=50):
        x_axis = [len(x) for x in self.readToContig.values()]
        plt.hist(x_axis, bins)
        plt.xlabel("Alignments per read")
        plt.ylabel("Number of reads")
        plt.savefig(title+".png")
        plt.clf()

