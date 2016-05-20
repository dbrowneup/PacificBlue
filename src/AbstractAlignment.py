#Re-written by Dan Browne April 2016

import sys

class AbstractAlignment:
    def __init__(self, aln, fileFormat="m4"):
        if fileFormat == "m5":
            self.qName = str(aln[0])
            self.qLength = int(aln[1])
            self.qStart = int(aln[2])
            self.qEnd = int(aln[3])
            self.qStrand = int(aln[4])
            self.tName = int(aln[5])
            self.tLength = int(aln[6])
            self.tStart = int(aln[7])
            self.tEnd = int(aln[8])
            self.tStrand = int(aln[9])
            self.score = int(aln[10])
            self.numMatch = int(aln[11])
            self.numMismatch = int(aln[12])
            self.numIns = int(aln[13])
            self.numDel = int(aln[14])
            self.mapQV = int(aln[15])
            self.qAlignedSeq = str(aln[16])
            self.matchPattern = str(aln[17])
            self.tAlignedSeq = str(aln[18])
        elif fileFormat == "m4":
            self.qName = str(aln[0])
            self.tName = int(aln[1])
            self.score = int(aln[2])
            self.percentSimilarity = float(aln[3])
            self.qStrand = int(aln[4])
            self.qStart = int(aln[5])
            self.qEnd = int(aln[6])
            self.qLength = int(aln[7])
            self.tStrand = int(aln[8])
            self.tStart = int(aln[9])
            self.tEnd = int(aln[10])
            self.tLength = int(aln[11])
            self.mapQV = int(aln[12])
        elif fileFormat == "m1":
            self.qName = str(aln[0])
            self.tName = int(aln[1])
            self.qStrand = int(aln[2])
            self.tStrand = int(aln[3])
            self.score = int(aln[4])
            self.percentSimilarity = float(aln[5])
            self.tStart = int(aln[6])
            self.tEnd = int(aln[7])
            self.tLength = int(aln[8])
            self.qStart = int(aln[9])
            self.qEnd = int(aln[10])
            self.qLength = int(aln[11])
        else:
            sys.exit("File format unrecognizable!")
    


