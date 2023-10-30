# parse hic pairs to removed duplication inter-ligation, intra-ligation, and correlated bam file.
import sys
import Utils as ut
import os
import gzip
import pysam


SAM_SEP = "\031"
SAM_INTER_SAM_SEP = "\031NEXT_SAM\031"


class PAIRTOOLSOBJ():
    def __init__(self, string):
        # print("String:", string)
        Cont = string.split("\t")
        self.rname = Cont[0]
        self.chr1 = Cont[1]
        self.pos1 = int(Cont[2])
        self.chr2 = Cont[3]
        self.pos2 = int(Cont[4])
        self.strand1 = Cont[5]
        self.strand2 = Cont[6]
        self.pairtoolstype = Cont[7]
        self.sam1 = Cont[8]
        self.sam2 = Cont[9]
        if len(Cont) > 9:
            self.addcol = Cont[10:]

class PAIRFRAG():
    def __init__(self, chr1='!', start1=0, end1=0, chr2='!', start2=0, end2=0, pet=0):
        self.frag = 0
        if chr1 != '!':
            self.frag = 1
        if chr2 != '!':
            self.frag = 2

        self.chr1 = chr1
        self.start1 = start1
        self.end1 = end1

        self.chr2 = chr2
        self.start2 = start2
        self.end2 = end2

        self.pet = pet

class DUPMachine():
    def __init__(self, name, chr_list):
        self.name = name
        self.duprec = dict()
        for chr in chr_list:
            self.duprec[chr] = set()
        self.dupcount = 0

    def _checkDup(self, new_frag):
        dupCount = 0
        if new_frag.chr1 in self.duprec.keys():
            if new_frag.start1 in self.duprec[new_frag.chr1]:
                dupCount = dupCount + 1
            else:
                self.duprec[new_frag.chr1].add(new_frag.start1)
        if new_frag.chr2 in self.duprec.keys():
            if new_frag.start2 in self.duprec[new_frag.chr2]:
                dupCount = dupCount + 1
            else:
                self.duprec[new_frag.chr2].add(new_frag.start2)

        if dupCount >= new_frag.frag:
            self.dupcount = self.dupcount+1
            return True
        else:
            return False

class HiCMachine():
    def __init__(self, file, chunksize=10000, prefix='', selfbp=0):
        self.file = file
        self.chunksize = chunksize
        self.chunk = []
        self.header = []
        self.complete = 2
        self.count = 0
        self.out_pair_bam = True
        self.out_5_pos = True
        self.prefix = prefix
        self.header_complete = False
        self.selfbp = int(selfbp)
        self.samheader = pysam.AlignmentHeader()
        self.chr_list = dict()


        self.barp = True
        self.barp_list = dict()

        self.mkdup = True
        self.duprecord = {}
        self.selfduprecordcount = 0
        self.pairduprecordcount = 0
        self.rmdup = False

        if self.rmdup:
            self.mkdup = True

        self.out_chunk_bam = []
        self.out_chunk_pair = []

        self.paircount = 0
        self.selfcount = 0
        self.nonecount = 0

    def _stringAttr(self, String):
        Cont = String.split("\t")
        # @SQ
        AttName = Cont[0]
        # SN:chr2L LN:12345
        AttDict = {}
        for Att in Cont[1:]:
            if Att.split(":")[0] == "LN":
                AttDict[Att.split(":")[0]] = int(Att.split(":")[1])
            else:
                AttDict[Att.split(":")[0]] = Att.split(":")[1]
        return AttName, AttDict
            
    def _samHeader(self):
        header = {}
        for Head in self.header:
            if Head.startswith("#samheader"):
                AttName, AttDict = self._stringAttr(Head.strip("#samheader: "))
                AttName = AttName.strip("@")
                if AttName in header.keys():
                    # print(type(AttName))
                    header[AttName].append(AttDict)
                else:
                    header[AttName] = [AttDict]
        return header

    def _getFPPos(self, sam):
        CigarList = sam.cigartuples
        FindCR = False

        # if read reverse, 5 primer would be the add all the base pairs after last reference consumes.
        if sam.is_reverse:
            if sam.reference_end != None:
                TotalLen = 0
                for i in range(len(CigarList)-1, 0, -1):
                    CigarTuple = CigarList[i]
                    if CigarTuple[0] not in ut.REFERENCECONSUMERPYSAMTUPLE:
                        TotalLen = TotalLen + CigarTuple[1]
                    else:
                        break
                return sam.reference_end + TotalLen

            else:
                return sam.reference_start
        else:
            if sam.reference_end != None:
                TotalLen = 0
                for i in range(len(CigarList)):
                    CigarTuple = CigarList[i]
                    if CigarTuple[0] not in ut.REFERENCECONSUMERPYSAMTUPLE:
                        TotalLen = TotalLen + CigarTuple[1]
                    else:
                        break
                return sam.reference_start - TotalLen
            else:
                return sam.reference_start
    
    def _samStrToAlign(self, samstr, ptype):
        AlignSegList = []
        FivePos = 0
        FiveChr = '!'
        if samstr != ".":
            for samrecord in samstr.split(SAM_INTER_SAM_SEP):
                AlignSeg = pysam.AlignedSegment().fromstring(samrecord.replace(SAM_SEP, "\t"), self.samheader)
                AlignSeg.set_tag(tag='LT', value=ptype, value_type='A')
                if not AlignSeg.is_secondary:
                    FivePos = self._getFPPos(AlignSeg)
                    FiveChr = AlignSeg.reference_name
                AlignSegList.append(AlignSeg)
        return AlignSegList, FiveChr, FivePos
    
    def _writeOutputFile(self, mode='w'):
        ## write header for pair
        self.outpair = gzip.open(self.prefix+"Pairs.gz", mode+'t')
        self.header.append('#samheader: @PG\tID:parsepairs\tPN:parsepairs\tVN:0.1.1\tCL:'+" ".join(sys.argv))
        self.header.append('#command: '+" ".join(sys.argv))

        for Line in self.header:
            if not Line.startswith("#columns"):
                self.outpair.write(Line+"\n")
            else:
                OutputHeaderTmp = Line.replace("#columns: ", "").split(" ")
                OutputHeaderTmp.remove('sam1')
                OutputHeaderTmp.remove('sam2')
                OutputHeader = OutputHeaderTmp[0:8]
                OutputHeader.append("ligatetype")
                if self.out_5_pos:
                    OutputHeader.append("fppos1")
                    OutputHeader.append("fppos2")
                if len(OutputHeaderTmp) > 8:
                    OutputHeader = OutputHeader + OutputHeaderTmp[8:]
                if self.barp:
                    OutputHeader.append("barpID")
                if self.mkdup and not self.rmdup:
                    OutputHeader.append("dupMark")
                self.outpair.write("#columns: "+" ".join(OutputHeader)+"\n")

            if Line.startswith("#chromsize"):
                Line = Line.replace("#chromsize: ", "").split(" ")
                self.chr_list[Line[0]] = int(Line[1])

        self.samheader = self.samheader.from_dict(self._samHeader())

        if self.out_pair_bam:
            self.outpb = pysam.AlignmentFile(self.prefix+"Pairs.bam", "wb", header=self._samHeader())
            # print(self.outpb.header)

    def _closeFiles(self):
        self.outpair.close()
        if self.out_pair_bam:
            self.outpb.close()

    def _readRecord(self):
        Count = len(self.chunk)
        while Count < self.chunksize:
            line = self.file.readline().strip()
            if not line:
                self.complete = 3
                break
            self.chunk.append(line)
            Count = Count + 1
        self.count = self.count + Count
        sys.stderr.write("Processing "+str(self.count)+"\n")
    
    def _outputChunk(self):
        self.outpair.writelines(self.out_chunk_pair)
        for i in self.out_chunk_bam:
            # print(i)
            self.outpb.write(i)
        self.out_chunk_bam = []
        self.out_chunk_pair = []

    def _parsePairRecord(self, Pair):
        PairObj = PAIRTOOLSOBJ(Pair)
        OutPairList = [PairObj.rname, PairObj.chr1, str(PairObj.pos1), PairObj.chr2, str(PairObj.pos2), PairObj.strand1, PairObj.strand2, PairObj.pairtoolstype]

        # paired interaction, inter ligation
        FivePos = ['0', '0']
        FiveChr = ['!', '!']
        SamPair = [PairObj.sam1, PairObj.sam2]
        RecordDup = False
        WriteRecord = True
        RecordType = 'N'

        AlignSegList = [[], []]
        if PairObj.pairtoolstype in {"UU", "UR", "RU"}:
            if PairObj.chr1 == PairObj.chr2 and (PairObj.pos2 - PairObj.pos1) < self.selfbp:
                self.selfcount = self.selfcount + 1
                OutPairList.append("self")
                RecordType='S'
                for i in range(2):
                    samstr = SamPair[i]
                    AlignSeg, SamFiveChr, SamFivePos = self._samStrToAlign(samstr=samstr, ptype='S')
                    FivePos[i] = str(SamFivePos)
                    FiveChr[i] = SamFiveChr
                    AlignSegList[i] = AlignSeg
            else:
                self.paircount = self.paircount + 1
                OutPairList.append("pair")
                RecordType='P'
                for i in range(2):
                    samstr = SamPair[i]
                    AlignSeg, SamFiveChr, SamFivePos = self._samStrToAlign(samstr=samstr, ptype='P')
                    FivePos[i] = str(SamFivePos)
                    FiveChr[i] = SamFiveChr
                    AlignSegList[i] = AlignSeg
        # for one end uniquely aligned, take Uniquely mapped end
        elif PairObj.pairtoolstype in {"UN", "UM", "RN", "RM"}:
            OutPairList.append("self")
            RecordType='S'
            self.selfcount = self.selfcount + 1
            samstr = SamPair[0]
            AlignSeg, SamFiveChr, SamFivePos = self._samStrToAlign(samstr=samstr, ptype='S')
            FivePos[0] = str(SamFivePos)
            FiveChr[0] = SamFiveChr
            AlignSegList[0] = AlignSeg
        elif PairObj.pairtoolstype in {"NU", "MU", "NR", "MR"}:
            OutPairList.append("self")
            RecordType='S'
            self.selfcount = self.selfcount + 1
            samstr = SamPair[1]
            AlignSeg, SamFiveChr, SamFivePos = self._samStrToAlign(samstr=samstr, ptype='S')
            FivePos[1] = str(SamFivePos)
            FiveChr[1] = SamFiveChr
            AlignSegList[1] = AlignSeg
        else:
            OutPairList.append("none")
            RecordType='N'
            self.nonecount = self.nonecount + 1
            for i in range(2):
                samstr = SamPair[i]
                AlignSeg, SamFiveChr, SamFivePos = self._samStrToAlign(samstr=samstr, ptype='N')
                FivePos[i] = str(SamFivePos)
                FiveChr[i] = SamFiveChr
                AlignSegList[i] = AlignSeg

        if self.out_5_pos:
            OutPairList = OutPairList + FivePos
        ## Add addition columns 
        OutPairList = OutPairList + PairObj.addcol
        if self.barp:
            DupName = ut.getBarpCode(PairObj.rname)
            OutPairList.append(DupName)
        else:
            DupName = 'All'
        
        if self.mkdup:
            if DupName in self.duprecord.keys():
                if self.duprecord[DupName]._checkDup(PAIRFRAG(chr1=FiveChr[0], start1=FivePos[0], chr2=FiveChr[1], start2=FivePos[1])):
                    if RecordType=='S':
                        self.selfduprecordcount = self.selfduprecordcount + 1
                    elif RecordType=='P':
                        self.pairduprecordcount = self.pairduprecordcount + 1
                    RecordDup = True
            else:
                self.duprecord[DupName] = DUPMachine(DupName, self.chr_list.keys())
                if self.duprecord[DupName]._checkDup(PAIRFRAG(chr1=FiveChr[0], start1=FivePos[0], chr2=FiveChr[1], start2=FivePos[1])):
                    if RecordType=='S':
                        self.selfduprecordcount = self.selfduprecordcount + 1
                    elif RecordType=='P':
                        self.pairduprecordcount = self.pairduprecordcount + 1
                    RecordDup = True
            if self.rmdup:
                if RecordDup:
                    WriteRecord = False
            else:
                if RecordDup:
                    OutPairList.append("dup")
                    if self.out_pair_bam:
                        for Aligns in AlignSegList:
                            if Aligns != []:
                                for Align in Aligns:
                                    Align.flag = Align.flag + 1024
                else:
                    OutPairList.append("nondup")


        if WriteRecord:
            self.out_chunk_pair.append("\t".join(OutPairList)+"\n")
            if self.out_pair_bam:
                for Aligns in AlignSegList:
                    if Aligns != []:
                        for Ali in Aligns:
                            self.out_chunk_bam.append(Ali)

    def _outStatistic(self):
        sys.stdout.write("\nProcessed: "+str(self.count)+" reads.\n")
        sys.stdout.write("Processed: "+str(self.paircount)+" pairs.\n")
        sys.stdout.write("Processed: "+str(self.selfcount)+" self-ligate fragments.\n")
        sys.stdout.write("Processed: "+str(self.nonecount)+" none reads.\n")
        sys.stdout.write("Processed: "+str(self.selfduprecordcount)+" duplicated self-ligate records.\n")
        sys.stdout.write("Processed: "+str(self.pairduprecordcount)+" duplicated pair records.\n")

    def work(self):
        while self.complete > 1:
            self._readRecord()
            if self.header_complete:
                while self.chunk:
                    Pair = self.chunk.pop(0)
                    if Pair.startswith("#"):
                        self.header.append(Pair.decode('utf-8'))
                    else:
                        self._parsePairRecord(Pair=Pair)
                self._outputChunk()
                if self.complete == 3:
                    break
            else:
                if self.complete == 3:
                    break
                while self.chunk:
                    Pair = self.chunk.pop(0)
                    if Pair.startswith("#"):
                        self.header.append(Pair)
                    else:
                        self.chunk.insert(0, Pair)
                        self._writeOutputFile()
                        self.header_complete=True
                        self.count = 0
                        break
        self._closeFiles()
        self._outStatistic()
                


def main():
    UsageList=[
        "Options:",
        "  -i [FILE]    Input pairs file produced by paritools parse.",
        "  -s [INT]     Self ligation base pair[3000].",
        "  -3 [INT]     Extend base pair from 3 prime end[0].", 
        "  -5 [INT]     IDB name, String before .idb as default[0].",
        "  -o [STR]     Output prefix."
        "  -c [INT]     Process chunk[10000]."
    ]

    ShortParaArgList=['i', 's', '3', '5', 'o', 'c']
    LongParaArgList=["input=", "selfbp=", "extbp3=", "extbp5=", "output=", "chunk="]
    ShortNonParaList=['v']
    LongNonParaList=["verbose"]
    Args = ut.parseOpt(sys.argv, ShortParaArgList=ShortParaArgList, 
                       LongParaArgList=LongParaArgList, 
                       UsageList=UsageList, 
                       LongNonParaList=LongNonParaList, 
                       ShortNonParaList=ShortNonParaList)
    
    FileInput, Selfbp, Extbp3, Extbp5, OutPrefix, Chunk = Args[0]

    if not os.path.exists(FileInput):
        sys.stderr.write("\nInput file " + FileInput + " not found...\n")
        ut.usagePrint(UsageList=UsageList)

    if Selfbp == '':
        Selfbp = 3000 
    if Extbp3 == '':
        Extbp3 = 0
    if Extbp5 == '':
        Extbp5 = 0
    if OutPrefix == '':
        OutPrefix = ""
    else:
        OutPrefix = OutPrefix + "."
    if Chunk == '':
        Chunk = 10000
    sys.stderr.write("\nInput:\t"+FileInput+
                     "\nExt 3 primer:\t"+Extbp3+" bp"+
                     "\nExt 3 primer:\t"+Extbp5+" bp"+
                     "\nSelf ligation:\t"+Selfbp+" bp"+
                     "\nOutput:\t"+OutPrefix+"All.pairs"+
                     "\n\n")
    
    if ut.isGzipFile(FileInput):
        Input = gzip.open(FileInput, 'rt')
    else:
        Input = open(FileInput, 'r', encoding='utf-8')
    
    WorkMachine = HiCMachine(Input, prefix=OutPrefix, selfbp=Selfbp)
    WorkMachine.work()
    with open(OutPrefix+"ParsePairs.stat", 'w', encoding='utf-8') as F:
        F.write("Read: "+str(WorkMachine.count)+"\n")
        F.write("Pair: "+str(WorkMachine.paircount)+"\n")
        F.write("Self: "+str(WorkMachine.selfcount)+"\n")
        F.write("None: "+str(WorkMachine.nonecount)+"\n")
        F.write("DuplicateSe: "+str(WorkMachine.selfduprecordcount)+"\n")
        F.write("DuplicatePe: "+str(WorkMachine.pairduprecordcount)+"\n")

if __name__ == "__main__":
    main()