#!/usr/bin/env python
#
# Create a list of commands for each shard
#
# ex) 1k bp block size for last run
# $ create_tasks.py -ff arctic_10.faa -if arctic_10.faa.idx -of arctic_10.faa.task -b 1000000 -dn arctic_10.faa.db
#
# Seung-Jin Sul (ssul@lbl.gov)
#
import argparse

# BLASTN_DEFAULTS = " -evalue 1e-30 -perc_identity 90 -word_size 45 -task megablast -outfmt 0 -show_gis -num_alignments 5 -num_descriptions 5 -dust yes -soft_masking true "
REFSEQ_MICROBIAL = "/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.microbial/refseq.microbial"
NT_BBDEDUPE_BBMASKED_FORMATTED = "/dna/shared/rqc/ref_databases/ncbi/CURRENT/nt/bbtools_dedupe_mask"
REFSEQ_MICROBIAL2 = "/scratch/rqc/refseq.microbial/refseq.microbial"
NT_BBDEDUPE_BBMASKED_FORMATTED2 = "/scratch/rqc/bbtools_dedupe_mask/nt_bbdedupe_bbmasked_formatted"

# example command list
# BLASTN_CMD = "module unload blast+/2.2.28; module load blast+/2.2.28; blastn "
BLASTN_CMD = " blastn "
# BLASTCMD = " %s -db %s -num_threads 8 -evalue 1e-30 -perc_identity 90 -word_size 45 -task megablast -outfmt 7 -show_gis -dust yes -soft_masking true -out %s "
BLASTCMD = " %s -db %s -num_threads 8 -evalue 1e-30 -perc_identity 90 -word_size 45 -task megablast -outfmt 7 -show_gis -dust yes -soft_masking true >> %s 2>> %s.log "
# blastOutfileName = "/global/projectb/scratch/sulsj/2014.06.18-blastmq-test/test_blastn_refseq_${USER}_${HOSTNAME}.out"
# blastnCmd = BLASTCMD % (BLASTN_CMD, REFSEQ_MICROBIAL2, blastOutfileName, blastOutfileName)
# blastOutfileName = "/global/projectb/scratch/sulsj/2014.06.18-blastmq-test/test_blastn_nt_${USER}_${HOSTNAME}.out"
# blastnCmd = BLASTCMD % (BLASTN_CMD, NT_BBDEDUPE_BBMASKED_FORMATTED2, blastOutfileName, blastOutfileName)

# TASKFORMAT = " module unload blast+/2.2.28; module load blast+/2.2.28; /global/projectb/scratch/sulsj/2014.06.18-blastmq-test/fasta_reader.py -i %(fastaFile)s -s %(startOffset)s -e %(endOffset)s | %(blastCmd)s ::0\n"
TASKFORMAT = " module unload blast+; module load blast+; /global/projectb/scratch/sulsj/2014.06.18-blastmq-test/fasta_reader.py -i %(fastaFile)s -s %(startOffset)s -e %(endOffset)s | %(blastCmd)s ::0\n"
# ./fasta_reader.py -i arctic_10.faa -s 0 -e 37375 | shifter --image=bryce911/lastal:869 lastal  -E0.00001 -f blasttab arctic_10.faa.db
LASTTASKFORMAT = "fasta_reader.py -i %(fastaFile)s -s %(startOffset)s -e %(endOffset)s | shifter --image=bryce911/lastal:869 lastal -E0.00001 -f blasttab %(db)s\n"
STARTENDFORMAT = "%(startOffset)s\t%(endOffset)s\n"

if __name__ == '__main__':
    ## usage = "python create_tasks.py -i indexFile -o outputTaskFileName -s blockSize"
    desc = "shard creator"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-if", "--index_file",
                        dest="indexFile",
                        help="Fasta index file",
                        required=True)
    parser.add_argument("-ff", "--fasta_file",
                        dest="fastaFile",
                        help="Fasta file",
                        required=True)
    parser.add_argument("-of", "--output_file",
                        dest="outTaskFileName",
                        help="Output task file name",
                        required=True)
    parser.add_argument("-bs", "--block_size",
                        dest="blockSize",
                        help="Query block size in byte",
                        required=True,
                        type=int)
    # parser.add_argument("-dn", "--ref_name",
    #                     dest="refName",
    #                     help="Ref database name: refseq or nt",
    #                     required=True,
    #                     type=str)
    options = parser.parse_args()

    with open(options.outTaskFileName, 'w') as TaskFH:
        with open(options.indexFile, 'r') as indexFH:
            ##
            ## index file line format: start_floc, seq_len, seq_num
            ##
            accSeqLen = 0
            blockStart = 0
            blockEnd = 0
            l = indexFH.readline()
            lastBlockEndRead = 0
            blockNum = 0

            while l:
                startLoc, endLoc, seqLen, seqNum = l.split()
                lastBlockEndRead = endLoc

                if accSeqLen == 0: blockStart = startLoc
                accSeqLen += int(seqLen)

                if accSeqLen >= options.blockSize:
                    blockEnd = endLoc

                    # print blockStart, blockEnd, accSeqLen, blockNum, l
                    # TaskFH.write(TASKFORMAT % {"fastaFile": options.fastaFile, "startOffset": blockStart, "endOffset": blockEnd, "blastCmd": blastnCmd})
                    # TaskFH.write(LASTTASKFORMAT % {"fastaFile": options.fastaFile,
                    #                                "startOffset": blockStart,
                    #                                "endOffset": blockEnd,
                    #                                "db": options.refName})
                    TaskFH.write(STARTENDFORMAT % {"startOffset": blockStart,
                                                   "endOffset": blockEnd})
                    accSeqLen = 0
                    blockNum += 1

                l = indexFH.readline()

            if int(blockEnd) < int(lastBlockEndRead):
                # print blockStart, lastBlockEndRead, accSeqLen, blockNum, l
                # TaskFH.write(TASKFORMAT % {"fastaFile": options.fastaFile, "startOffset": blockStart, "endOffset": lastBlockEndRead, "blastCmd": blastnCmd})
                TaskFH.write(
                    # LASTTASKFORMAT % {"fastaFile": options.fastaFile,
                    #                   "startOffset": blockStart,
                    #                   "endOffset": lastBlockEndRead,
                    #                   "db": options.refName})
                    STARTENDFORMAT % {"startOffset": blockStart,
                                      "endOffset": blockEnd})
                blockNum += 1

## EOF
