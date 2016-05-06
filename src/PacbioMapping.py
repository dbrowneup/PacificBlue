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
        self.alignments = set([AbstractAlignment(x.split(' '), fileFormat) for x in self.alignments]) # parses alignments into AbstractAlignment objects
        self.readToContig = {x.qName: set([]) for x in self.alignments if x.qName not in self.readToContig} # parses qNames into dictionary
        #parse alignments into dictionary
        for x in self.alignments:
            self.readToContig[x.qName].add(x)
        print "Number of mapping reads:", len(self.readToContig)
        print "Number of alignments:", len(self.alignments)
        #Filter out reads with only 1 alignment
        self.filter_reads()
        print "Finished loading PacBio mapping:", str(datetime.now())

    def filter_reads(self):
        filtered_reads = set([])
        for k,v in self.readToContig.items():
            if len(v) == 1:
                filtered_reads.add(k)
                del self.readToContig[k]
        self.alignments = set([x for x in self.alignments if x.qName not in filtered_reads])
        print "Number of mapping reads filtered out:", len(filtered_reads)
        print "Number of mapping reads remaining:", len(self.readToContig)

    def read_mapping_frequency(self, title, bins=50):
        x_axis = [len(x) for x in self.readToContig.values()]
        plt.hist(x_axis, bins)
        plt.xlabel("Alignments per read")
        plt.ylabel("Number of reads")
        plt.savefig(title+".png")
        plt.clf()

