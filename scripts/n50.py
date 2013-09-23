#!/usr/bin/env python

# calculate N50 from assembled contig fasta file

# If you have R, this script will generate a PDF histogram of the contig lengths.


import commands
import sys
import os
from Bio import SeqIO

lengths = []
for record in SeqIO.parse(open(sys.argv[1]),'fasta'):
	lengths.append(len(record.seq))

lengths.sort(reverse=True)
teoN50 = sum(lengths)/2.0
testSum = 0
N50 = 0
for con in lengths:
	testSum += con
	if teoN50 < testSum:
		N50 = con
		break
print "# of contigs:", len(lengths)
print "N50 = %s" % N50

if os.path.exists(commands.getoutput('which r')):
    
    with open('all_lengths.txt', 'w') as handle:
      handle.write('\n'.join(str(i) for i in lengths))

    histogram_calculator = """
    data <- read.csv('all_lengths.txt', header=F)
    hist(as.matrix(data), main='N50 = %s', xlab='length', ylab='frequency')
    abline(v=%s, col='red')
    """ % (n50, n50)

    with open('histogram.r', 'w') as handle:
      handle.write(histogram_calculator)
    
    commands.getoutput('r --slave < histogram.r')
    commands.getoutput('rm histogram.r')
    
    if len(sys.argv) == 3:
        commands.getoutput(' '.join(['mv', 'Rplots.pdf', sys.argv[2]]))
