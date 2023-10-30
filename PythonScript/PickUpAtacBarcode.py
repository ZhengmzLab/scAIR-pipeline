# This script correlate two identifier based on the barcode order.
import sys
import Utils as ut




def main():
    UsageList=[
        "Options:",
        "  -a FILE    First IDB file.",
        "  -b FILE    Second IDB file.", 
        "  -n STR     IDB name, String before .idb as default.",
        "  -o FILE    Output prefix."
    ]

    ShortParaArgList=['a', 'b', 'n', 'o']
    LongParaArgList=["reference=", "query=", "idb=", "output="]
    ShortNonParaList=['v']
    LongNonParaList=["verbose"]
    Args = ut.parseOpt(sys.argv, ShortParaArgList=ShortParaArgList, 
                       LongParaArgList=LongParaArgList, 
                       UsageList=UsageList, 
                       LongNonParaList=LongNonParaList, 
                       ShortNonParaList=ShortNonParaList)
    Fidb, Sidb, IdbName, OutPrefix = Args[0]
    if Fidb.strip('.idb').split(".")[-1] != Sidb.strip('.idb').split(".")[-1]:
        sys.stderr.write("You are trying to correlating two different type of identifier...\n")
    
    if IdbName == "":
        IdbName = Fidb.strip('.idb').split(".")[-1]

    if OutPrefix == "":
        OutPrefix = ".".join(Fidb.strip('.idb').split(".")[:])+"."+".".join(Sidb.strip('.idb').split(".")[:])+".correlate.tsv"
    else:
        OutPrefix = OutPrefix + ".correlate.tsv"

    FidbFile = open(Fidb, "r", encoding='utf-8')
    SidbFile = open(Sidb, "r", encoding='utf-8')

    FirstBarList = {}
    while True:
        line = FidbFile.readline().strip().split()
        if not line:
            break
        if not FirstBarList.get(line[0]):
            FirstBarList[line[0]] = line[1]
        else:
            sys.stderr.write("One barcode corresponding two identifier.\n")
            FirstBarList[line[0]] = line[1]

    SecondBarList = {}
    while True:
        line = SidbFile.readline().strip().split()
        if not line:
            break
        if not SecondBarList.get(line[0]):
            SecondBarList[line[0]] = line[1]
        else:
            sys.stderr.write("One barcode corresponding two identifier.\n")
            SecondBarList[line[0]] = line[1]

    AllKeys = set(list(FirstBarList.keys()) + list(SecondBarList.keys()))
    with open(OutPrefix, "w", encoding='utf-8') as F:
        F.write("#"+".".join(Fidb.strip('.idb').split(".")[:])+"\t"+".".join(Sidb.strip('.idb').split(".")[:])+"\n")
        for key in AllKeys:
            FirstColumn="-"
            SecondColumn="-"
            if FirstBarList.get(key):
                FirstColumn=IdbName+"_"+FirstBarList[key]
            if SecondBarList.get(key):
                SecondColumn=IdbName+"_"+SecondBarList[key]
            F.write(FirstColumn + "\t" + SecondColumn + "\n")


if __name__ == "__main__":
    main()



