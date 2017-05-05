#Re-written by Dan Browne on 04/27/16


import numpy as np
import itertools as it

class PacbioSubgraph():

    def __init__(self, alignments, cov_cutoff, fraction):
        self.mapping = alignments
        #Create array of read aligments by read length
        self.readArray = self.read_stack()
        #Mark coordinates of alignments on read
        for i, a in enumerate(self.mapping):
            self.mark_coords(i, a.qStart, a.qEnd)
        #Calculate per-base coverage of read
        self.covArray = self.readArray.sum(axis=0)
        #Filter out alignments in repetitive regions
        self.mapping = set([a for a in self.mapping if self.uniquity_score(a) > fraction])
    
    def __enter__(self):
        return self

    def __exit__(self, exctype, value, traceback):
        if exctype is None:
            return

    def make_connects(self):
        #Find 5' and 3' overhangs
        overhang5, overhang3 = self.find_overhangs()
        #Connect 5' and 3' overhangs
        connects = list()
        possible_connects = list(it.product(overhang3, overhang5))
        connects = it.imap(self.test_connect, possible_connects)
        return [c for c in connects if c is not None]
    
    def read_stack(self):
        stack_size = len(self.mapping)
        read_size = next(iter(self.mapping)).qLength
        stack = np.zeros((stack_size, read_size))
        return stack
    
    def mark_coords(self, n, qStart, qEnd):
        for i in xrange(qStart, qEnd):
            self.readArray[n][i] = 1

    def uniquity_score(self, a):
        S = sum((float(1)/self.covArray[i] for i in xrange(a.qStart, a.qEnd)))
        L = float(a.qEnd - a.qStart)
        MAPU = S/L
        return MAPU
    
#    def find_repeats(self, cutoff):
#        repeat_coords = []
#        for i in range(len(self.covArray)):
#            if self.covArray[i] > cutoff:
#                repeat_coords.append(i)
#        return repeat_coords

#    def filter_repeat(self, a, fraction):
#        align_coords = set(range(a.qStart, a.qEnd))
#        bp_in_repeat = align_coords.intersection(self.repeat_coords)
#        if len(bp_in_repeat) >= (fraction * len(align_coords)):
#            return True

    def find_overhangs(self):
        overhang5 = list()
        overhang3 = list()
        for a in self.mapping:
            if a.qStart > a.tStart and (a.qLength - a.qEnd) < (a.tLength - a.tEnd):
#                print "5' overhang"
                overhang5.append(a)
            elif a.qStart < a.tStart and (a.qLength - a.qEnd) > (a.tLength - a.tEnd):
#                print "3' overhang"
                overhand3.append(a)
            elif a.qStart >= a.tStart and (a.qLength - a.qEnd) >= (a.tLength - a.tEnd):
#                print "Double overhang or blunt"
                return
            elif a.qStart <= a.tStart and (a.qLength - a.qEnd) <= (a.tLength - a.tEnd):
#                print "No overhang or blunt"
                return
            else:
                print "Unknown overhang!"
        return (overhang5, overhand3)
    
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

