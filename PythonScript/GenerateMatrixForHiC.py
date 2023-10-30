import sys
import Utils as ut
import os


def main():
    UsageList=[
        "Options:",
        "  -i [File]   Input file [contain barp ID and fragment count column]."
        "  -b [INT]    BARP ID column [1-base, default: 4].",
        "  -c [INT]    Count columnn [1-base, default: 5].", 
        "  -o [STR]    Output prefix."
    ]

    ShortParaArgList=['i', 'b', 'c', 'o']
    LongParaArgList=["input=", "barpidc=", "countc=", "output="]
    ShortNonParaList=['v']
    LongNonParaList=["verbose"]

    Args = ut.parseOpt(sys.argv, ShortParaArgList=ShortParaArgList, 
                    LongParaArgList=LongParaArgList, 
                    UsageList=UsageList, 
                    LongNonParaList=LongNonParaList, 
                    ShortNonParaList=ShortNonParaList)

    InputFile, BarpIDCol, CountCol, OutputPrefix = Args[0]
    Verbose = Args[1]

    if not os.path.exists(InputFile):
        sys.stderr.write("\nInput file " + InputFile + " not found...\n")
        ut.usagePrint(UsageList=UsageList)

    if not OutputPrefix == '':
        OutputPrefix = OutputPrefix + "."
    if not BarpIDCol == '':
        BarpIDCol = int(BarpIDCol)-1
    else:
        BarpIDCol = 3
    
    if not CountCol == '':
        CountCol = int(CountCol)-1
    else:
        CountCol = 4

    # open files
    FI = open(InputFile, "r")
    OFeature = open(OutputPrefix + "features.tsv", "w")
    OBarp = open(OutputPrefix + "barcodes.tsv", "w")
    OMatrix = open(OutputPrefix + "matrix.mtx", "w")
    
    MatrixMarketComment = "%%MatrixMarket matrix coordinate integer general\n" + "%GeneMatrix\tCL:"+" ".join(sys.argv) + "\n"
    OMatrix.write(MatrixMarketComment)

    FeatureList = {}
    CellList = {}

    FeatureCount = 1
    CellCount = 1
    LineCount = 0
    MatrixList = []

    while True:
        line = FI.readline().strip()
        if not line:
            break
        line = line.split()
        CurFeature = "-".join(line[:BarpIDCol])
        CurCell = line[BarpIDCol]
        LineCount = LineCount + 1

        if not FeatureList.get(CurFeature):
            FeatureList.update({CurFeature: FeatureCount})
            FeatureLine = CurFeature.replace("-", "\t")+"\n"
            OFeature.write(FeatureLine)
            FeatureIndex = FeatureCount
            FeatureCount = FeatureCount + 1
        else:
            FeatureIndex = FeatureList.get(CurFeature)

        if not CellList.get(CurCell):
            CellList.update({CurCell: CellCount})
            CellLine = CurCell+"\n"
            OBarp.write(CellLine)
            CellIndex = CellCount
            CellCount = CellCount + 1
        else:
            CellIndex = CellList.get(CurCell)
    
        FragCount = line[CountCol]
        MatrixLine = str(FeatureIndex)+"\t" + str(CellIndex)+"\t"+FragCount+"\n"
        MatrixList.append(MatrixLine)

        if (LineCount % 10000 == 0) and Verbose:
            sys.stderr.write("[processed] " + str(LineCount) + "\n")
        
    DimensionLine = str(FeatureCount-1)+"\t"+str(CellCount-1)+"\t"+str(LineCount)+"\n"
    OMatrix.write(DimensionLine)

    for FeatureCellCount in MatrixList:
        OMatrix.write(FeatureCellCount)
    
    OBarp.close()
    OFeature.close()
    OMatrix.close()
    FI.close()
    print("Generate Matrix Finished. Processed: " + str(LineCount) + ".")

if __name__ == "__main__":
    main()