#!/usr/bin/env python
# indexes fastq reads for shard_reader.py
# Jeff Froula <jlfroula@lbl.gov>
# 08/21/2019
#
import sys
import optparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

uid=0
recordnum=0
seqLen=0
start=0

usage = "fasta_indexer.py -i <fastq> -o <outIndexFile>"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-i", "--input", dest="infilename",
                      action="store", type="string", help="input fastq file", default=False)
parser.add_option("-o", "--output", dest="outfilename",
                      action="store", type="string", help="output index file")
(options, args) = parser.parse_args()

if options.infilename:
	inFileName = options.infilename
else:
	parser.error("Please set the input file name.")

f = open(options.outfilename, "w")

with open(options.infilename, "r") as handle:
    for title, seq, qual in FastqGeneralIterator(handle) :
        id_len = len(title)+1
        seq_len = len(seq)
        qual_len = len(qual)

        if seq_len != qual_len:
            print("sequence length not equal to quality score length for read {}".format(title))
            print("Please run reformat.sh from bbmap package to fix your fastq file")
            sys.exit(1)

        # "start" index starts from 0 and not 1. I also need 
        # to count the + line and each \n
        end = start + id_len + seq_len + qual_len + 4

        f.write("{}\t{}\t{}\t{}\n".format(start,end,seq_len,uid))

        start=end+1
        uid+=1
f.close()
