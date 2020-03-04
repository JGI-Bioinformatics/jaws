#!/usr/bin/env python
#
# Create a index file for input fasta file using bp block size
#
# ex)
# $ fasta_indexer.py -i arctic_10.faa -o arctic_10.faa.idx -d arctic_10.faa.def -s 0 -b 1
#
# Seung-Jin Sul (ssul@lbl.gov)
#
import optparse
import gzip
import os
from subprocess import Popen, call, PIPE


def openCompressed(filename, mode, compressFormat=None, **kw):
    """Open a filename which can be either compressed or plain file.
    @param compressFormat if None, an attempt will be made to autodetect format
    (currently by extension, only '.gz' and '.bz2' are recognized); if "none" - open as plain
    file, if "gzip" - open as gzip file."""
    cform = compressFormat
    if cform is None:
        cform = "none"
        if filename.endswith('.gz'):
            cform = "gzip"
        elif filename.endswith('.bz2'):
            cform = "bz2"
    # print("DEBUG: openCompressed(%s,%s,%s)" % (filename, mode, cform))
    k = kw.copy()
    if cform == "gzip":
        if "buffering" in k:
            k["bufsize"] = k["buffering"]
            del k["buffering"]
        return openGzip(filename, mode, **kw)
    elif cform == "bz2":
        k.setdefault("buffering",2**20)
        return bz2.BZ2File(filename, mode, **kw)
    elif cform == "none":
        k.setdefault("buffering",2**20)
        return open(filename, mode, **kw)
    else:
        raise ValueError(compressFormat)


def openGzip(filename, mode, compresslevel=6):
    compresslevel = int(compresslevel)

    if mode in ("w", "wb", "a", "ab"):
        if mode in ("w", "wb"):
            redir = ">"
        elif mode in ("a", "ab"):
            redir = ">>"
        # p = Popen("gzip -%s %s %s" % (compresslevel,redir,filename), shell=True, env=os.environ, bufsize=2**16, stdin=PIPE, close_fds=True)
        stdout = open(filename, mode)
        p = Popen(["gzip", "-%s" % compresslevel, "-c"], shell=False, env=os.environ,
                  bufsize=2**16, stdin=PIPE, stdout=stdout, close_fds=True)
        stdout.close()
        return PopenStdinProxy(p)

    elif mode in ("r", "rb"):
        cmd = ["gzip", "-cd", filename]
        if not os.path.isfile(filename):
            raise OSError("Input file does not exist: %s", filename)
        p = Popen(cmd, env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True)
        if p.returncode:
            raise CalledProcessError(str(cmd), p.returncode)
        return p.stdout

    else:
        raise ValueError("'openGzip()' - Unsupported 'mode': " + mode)


class FastaReader(object):
    """Class that supports an input iterator protocol for a FASTA file.
    Example that prints an exact copy of the input file:
    for rec in FastaReader(open('seq.fsa','r')).records():
        print rec.header(),
        for line in rec.seqLines():
            print line,
    Instead of rec.seqLines(), you can use the methods which post-process the
    raw sequences lines: seqChunks(), seqArrays(), sequence().
    """

    def __init__(self, infile):
        # if infile.endswith("gz"):
        #     infile = gzip.open(infile, 'r')
        # else:
        #     infile = open(infile, 'r')
        if not hasattr(infile, "readline"):
            infile = openCompressed(infile, 'r')
            # infile = gzip.open(infile, 'r', encoding='utf-8', errors='ignore')
            # infile = codecs.open(infile, 'r', encoding='utf-8', errors='ignore')
            # fd = gzip.open(infile, 'rb')
            # infile = io.BufferedReader(fd)
            self.ownInfile = True
        else:
            self.ownInfile = False

        print(self.ownInfile)

        self.infile = infile
        self.freshHdr = False
        self.maxLineLen = 0

    def records(self):
        infile = self.infile
        while True:
            if self.freshHdr:
                self.freshHdr = False
                yield self
                continue
            line = infile.readline()
            if type(line).__name__ == 'bytes':
                line = line.decode()

            if not line:
                return
            # skip blank lines
            elif line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                yield self

    def header(self):
        assert self.hdr.startswith('>')
        return self.hdr

    def getNCBI_Id(self):
        """Assume that header starts with '>gi|1234567|' and return the string id from second field."""
        return self.hdr.split('|', 2)[1]

    def getNCBI_GI(self):
        return int(self.getNCBI_Id())

    def getSimpleId(self):
        """Assume that header starts with '>string_no_spaces ' and return that string."""
        return self.hdr.split(None, 1)[0][1:]

    def seqLines(self):
        infile = self.infile
        while True:
            line = infile.readline()
            if type(line).__name__ == 'bytes':
                line = line.decode()
            if not line:
                break
            elif line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                self.freshHdr = True
                return
            self.maxLineLen = max(self.maxLineLen, len(line) - 1)
            yield line

    def seqChunks(self, chunkSize):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
            if seq.tell() >= chunkSize:
                yield seq.getvalue()
                seq.close()
                seq = StringIO()
        if seq.tell() > 0:
            yield seq.getvalue()
        seq.close()

    # def seqArrays(self,chunkSize):
    # for s in self.seqChunks(chunkSize):
    # yield numpy.fromstring(s,dtype='S1')

    def sequence(self, format='str'):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
        s = seq.getvalue()
        seq.close()
        # if format == 'array':
        # s = numpy.fromstring(s,dtype='S1')
        return s

    def seqLen(self):
        n = 0
        for line in self.seqLines():
            n += len(line) - 1
            if not line.endswith("\n"):
                n += 1
        return n

    def lineLen(self):
        return self.maxLineLen

    def close(self):
        if self.ownInfile:
            self.infile.close()


version = '2.0'
verbose = False
inFileName = ''
outFileName = ''
defFileName = ''
deflineOption = 0
uidOption = 0
seqUID = 0

if __name__ == '__main__':
    # usage = "python fasta_indexer.py -i inFile -o outIndexFile -d outDeflineFile -s startNum -b full_or_part_defline"
    usage = "python fasta_indexer.py -i inFile -o outIndexFile"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-i", "--input", dest="infilename",
                      action="store", type="string", help="input fasta file", default=False)
    parser.add_option("-o", "--output", dest="outfilename",
                      action="store", type="string", help="output index file")
    # parser.add_option("-d", "--deflinefilename", dest="deflinefilename",
    #                   action="store", type="string", help="defline file name")
    # parser.add_option("-s", "--startno", dest="startno",
    #                   action="store", type="int", help="uid start number", default=1)
    # parser.add_option("-b", "--deflineopt", dest="deflineopt",
    #                   action="store", type="int", help="defline saving option: 0=part, 1=full")
    (options, args) = parser.parse_args()

    if options.infilename:
        inFileName = options.infilename
    else:
        parser.error("Please set the input file name.")

    if options.outfilename:
        outFileName = options.outfilename
    else:
        parser.error("Please set the output index file name.")

    # if options.startno is not None:
    #     seqUID = options.startno
    #
    # if options.deflinefilename and options.deflineopt is not None:
    #     defFileName = options.deflinefilename
    #     deflineOption = options.deflineopt
    # else:
    #     parser.error("Please set the output def file name.")

    numSeq = 0
    currLoc = 0
    outFile = open(outFileName, "w")
    # defFile = open(defFileName, "w")

    # for rec in FastaReader(fd).records():
    for rec in FastaReader(inFileName).records():
        defline = rec.header()

        ###
        ### Save start location of each query def line and the length
        ###
        loc = currLoc
        currLoc += len(defline)
        seqLen = 0
        for line in rec.seqLines():
            seqLen += len(line.rstrip("\n"))
            currLoc += len(line)
        numSeq += 1
        uid = seqUID

        outFile.write(str(loc) + "\t" + str(currLoc - 1) + "\t" + str(seqLen) + "\t" + str(uid) + "\n")

        # deflines = ''
        # if deflineOption == 0:  ## save part of defline before the first blank
        #     defline2 = defline.rstrip().split(" ")[0]
        # else:  ## save full defline
        #     defline2 = defline.rstrip()
        # defFile.write(str(uid) + "\t" + defline2 + "\n")

        seqUID += 1

    outFile.close()
    # defFile.close()
    print(" %s: total number of sequences = %s" % (outFileName, numSeq))

### EOF
