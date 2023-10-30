import sys
import getopt
import time
import gzip


# correspond to pysam tuple
REFERENCECONSUMERPYSAMTUPLE={0, 2, 3, 7, 8}
QUERYCONSUMERPYSAMTUPLE={0, 1, 4, 7, 8}


class IDB():
    def __init__(self, barorder=-1, idenorder=0, readcount=0, idenstr=""):
        self.barorder = barorder
        self.idenorder = idenorder
        self.readcount = readcount
        self.idenstr = idenstr

def isGzipFile(Input):
    try:
        with gzip.open(Input) as F:
            F.read(1)
        return True
    except IOError:
        return False

def usagePrint(UsageList):
    sys.stderr.write("\n")
    for i in range(len(UsageList)):
        sys.stderr.write(UsageList[i]+"\n")
    sys.exit(2)

def parseOpt(ArgList, ShortNonParaList=[], ShortParaArgList=[], LongNonParaList=[], LongParaArgList=[], UsageList=""):
    if len(ShortNonParaList)+len(ShortParaArgList) != len(LongNonParaList)+len(LongParaArgList):
        raise RuntimeError("Number of short arguments do not match long arguments.")
    try:
        opts, args = getopt.getopt(ArgList[1:], ":".join(ShortParaArgList)+":"+"".join(ShortNonParaList)+"h", LongParaArgList+LongNonParaList+["--help"])
    except getopt.GetoptError:
        usagePrint(UsageList=UsageList)
    RetArgList = [""]*len(ShortParaArgList)
    for i in range(len(ShortParaArgList)):
        for j in range(len(opts)):
            opt = opts[j][0]
            arg = opts[j][1]
            if opt.strip("-") in [ShortParaArgList[i], LongParaArgList[i].strip("=")]:
                if RetArgList[i] != "":
                    sys.stderr.write("Argument " + opt + "appeared twice, the later one will coverage the before one.\n")
                RetArgList[i] = arg

    RetNonArgList = [False]*len(ShortNonParaList)
    for i in range(len(ShortNonParaList)):
        for j in range(len(opts)):
            opt = opts[j][0]
            arg = opts[j][1]
            if opt.strip("-") in [ShortParaArgList[i], LongParaArgList[i]]:
                RetNonArgList[i] = True
    
    return (RetArgList, RetNonArgList)

def getBarpCode(Input):
    return Input.strip().split("|||")[2]

class FASTQ():
    def __init__(self, name, seq, comment, quality):
        self.name = name
        self.seq = seq
        self.comment = comment
        self.quality = quality
    def _output(self):
        return self.name+"\n"+self.seq+"\n"+self.comment+"\n"+self.quality+"\n"
    def _complete(self):
        if self.name==0:
            return 0
        else:
            return 1
        
