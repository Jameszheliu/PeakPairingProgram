from collections import defaultdict
from bisect import bisect_left
from datetime import datetime

import sys
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-l', action='store', type='int', dest='length',
                      help='Limitation between Left and Right Peaks in bp')
(options, args) = parser.parse_args()

if not args:
    print "********************************************************************************************************"
    print "This program takes a ChIP-Exo Peak file and does peak pairing"
    print "Example: python PeakPairing.py Initial -l 20"
    print "Left, Right peaks should be stored in .peak format as Initial_Left.Peak and Initial_Rifht.Peak"
    print "Peak format: chr   strand(+/-)   Start   End   Score"
    print "-l: The distance allowed between the left and right peak middle point"
    print "output file is stored in finalpeakpairs.txt"
    print "********************************************************************************************************"
    sys.exit(1)

if not options.length:
    print "Read Length is missing"
    sys.exit(1)
    
PeakLen = int(options.length)
Initial = args[0]

def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    after = myList[pos]
    return after
  
print "Exo Peak Pairing"
print "start time", datetime.now()

file1 = Initial + '_Right.Peak'

fh = open(file1, 'r')
Right = list()
fh.readline()
for line in fh:
    nline = line.strip().split('\t')
    for item in nline:
        Right.append(item)
fh.close()

num = len(Right)/5
Total = defaultdict(list)
Total2 = defaultdict(list)

#make dictionary
for i in range(0, num):
    Total[Right[i*5]].append(int((int(Right[i*5 + 2]) + int(Right[i*5 + 3]))/2))
    Total2[Right[i*5]].append(Right[i*5 + 4])

file2 = Initial + '_Left.Peak'

fh2 = open(file2, 'r')
fh2.readline()

oh = open('finalpeakpairs.txt', 'w')

for line in fh2:
    nline = line.strip().split('\t')
    myKey = nline[0]
    myNumber = int((int(nline[2]) + int(nline[3]))/2)
    
    if Total.has_key(myKey):
        myList = Total[myKey]
        closest = takeClosest(myList, myNumber)
            
        if (closest - myNumber) < PeakLen:
            idx = myList.index(closest)
            RightScore = Total2[myKey][idx]
            rline = '\t'.join([nline[0], str(myNumber-1), str(closest), nline[4], RightScore + '\n'])
            oh.write(rline)

oh.close()
fh2.close()

print "Done with PeakPairing", datetime.now()
            
