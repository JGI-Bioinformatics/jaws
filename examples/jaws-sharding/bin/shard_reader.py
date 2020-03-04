#!/usr/bin/env python
#
# Catting a fasta shard to stdout
#
# ex)
# $ fasta_reader.py -i arctic_10.faa -s 0 -e 1458736
#
# Seung-Jin Sul (ssul@lbl.gov)
#
import argparse

if __name__ == '__main__':
    ## usage = "python fasta_reader.py -i inFile -s startOffset -e endOffset"
    desc = "shard reader"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-i", "--fasta", dest="fastaFile", help="Fasta query file", required=True)
    parser.add_argument("-s", "--start", dest="startOffset", help="start offset", required=True,
                        type=int)
    parser.add_argument("-e", "--end", dest="endOffset", help="end offset", required=True, type=int)
    options = parser.parse_args()
    with open(options.fastaFile, 'r') as FH:
        FH.seek(options.startOffset, 0)
        lastPos = 0
        while lastPos <= options.endOffset:
            print(FH.readline().rstrip())
            lastPos = FH.tell()

## EOF
