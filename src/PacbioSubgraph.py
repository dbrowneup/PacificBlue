#Re-written by Dan Browne on 04/27/16


import numpy as np
import itertools

class PacbioSubgraph():

    def __init__(self, read, alignments, cov_cutoff=3, fraction=0.5):
        self.mapping = alignments
        #Create array of read aligments by read length
        self.readArray = self.read_stack()
        #Mark coordinates of alignments on read
        stack_index = 0
        for a in self.mapping:
            self.mark_coords(stack_index, a.qStart, a.qEnd)
            stack_index += 1
        #Calculate per-base coverage of read
        self.covArray = self.readArray.sum(axis=0)
        #Filter out alignments in repetitive regions
        self.repeat_coords = set(self.find_repeats(cov_cutoff))
        self.mapping = set([a for a in self.mapping if not self.filter_repeat(a, fraction)])
        #Find 5' and 3' overhangs
        self.FivePOvhg = []
        self.ThrePOvhg = []
        for a in self.mapping:
            self.find_overhangs(a)
        #Connect 5' and 3' overhangs
        self.Connects = []
        possible_connects = list(itertools.product(self.ThrePOvhg, self.FivePOvhg))
        self.Connects = map(self.test_connect, possible_connects)
        self.Connects = [c for c in self.Connects if c is not None]
    
    def read_stack(self):
        stack_size = len(self.mapping)
        read_size = next(iter(self.mapping)).qLength
        stack = np.zeros((stack_size, read_size))
        return stack
    
    def mark_coords(self, n, qStart, qEnd):
        for i in range(qStart, qEnd):
            self.readArray[n][i] = 1
    
    def find_repeats(self, cutoff):
        repeat_coords = []
        for i in range(len(self.covArray)):
            if self.covArray[i] > cutoff:
                repeat_coords.append(i)
        return repeat_coords

    def filter_repeat(self, a, fraction):
        align_coords = set(range(a.qStart, a.qEnd))
        bp_in_repeat = align_coords.intersection(self.repeat_coords)
        if len(bp_in_repeat) >= (fraction * len(align_coords)):
            return True

    def find_overhangs(self, a):
        if a.qStart > a.tStart and (a.qLength - a.qEnd) < (a.tLength - a.tEnd):
#            print "5' overhang"
            self.FivePOvhg.append(a)
        elif a.qStart < a.tStart and (a.qLength - a.qEnd) > (a.tLength - a.tEnd):
#            print "3' overhang"
            self.ThrePOvhg.append(a)
        elif a.qStart >= a.tStart and (a.qLength - a.qEnd) >= (a.tLength - a.tEnd):
#            print "Double overhang or blunt"
            return
        elif a.qStart <= a.tStart and (a.qLength - a.qEnd) <= (a.tLength - a.tEnd):
#            print "No overhang or blunt"
            return
        else:
            print "Unknown overhang!"
    
    def test_connect(self, p):
        tp, fp = p
        if tp.qEnd <= fp.qStart:
            if (fp.qEnd - tp.qStart) > (tp.tLength - tp.tEnd + fp.tStart) <= (fp.qStart - tp.qEnd):
#                print "Read spans gap between scaffolds or connects blunt ends"
                gap_estimate = (fp.qStart - tp.qEnd) - (tp.tLength - tp.tEnd + fp.tStart)
                return (tp, fp, gap_estimate)
            elif (fp.qEnd - tp.qStart) > (tp.tLength - tp.tEnd + fp.tStart) > (fp.qStart - tp.qEnd):
#                print "Gap between alignments, possible overlap of scaffold ends"
                overlap_estimate = (fp.qStart - tp.qEnd) - (tp.tLength - tp.tEnd + fp.tStart)
                return (tp, fp, overlap_estimate)
            elif (fp.qEnd - tp.qStart) <= (tp.tLength - tp.tEnd + fp.tStart) > (fp.qStart - tp.qEnd):
#                print "Scaffold ends extend past alignment boundaries"
                return
            else:
                print "Unknown alignment condition 1!"
        elif fp.qEnd >= tp.qEnd > fp.qStart >= tp.qStart:
            if fp.tEnd < (fp.qEnd - tp.qStart) > (fp.tLength - fp.tStart):
#                print "Read spans overlap of scaffold ends"
                overlap_estimate = (fp.qEnd - tp.qStart) - (fp.tLength - fp.tStart + fp.tEnd)
                return (tp, fp, overlap_estimate)
            else:
#                print "Scaffolds extend beyond alignment boundaries"
                return
        elif tp.qEnd >= fp.qEnd > tp.qStart >= fp.qStart:
#            print "Inversion of expected alignment orientation"
            return
        elif fp.qEnd >= tp.qEnd > tp.qStart >= fp.qStart:
#            print "A1 contained within A2"
            return
        elif tp.qEnd >= fp.qEnd > fp.qStart >= tp.qStart:
#            print "A2 contained within A1"
            return
        elif tp.qStart >= fp.qEnd:
#            print "Gap between inverted alignments, unlikely scaffold overlap"
            return
        else:
#            return
            print "Unknown alignment condition 2!"

