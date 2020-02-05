from sys import argv
import re
import statistics

inputFile = argv[1]
outputFile = argv[2]

with open(inputFile) as chFile, open(outputFile, "w") as mscFile:
    header = chFile.readline()
    headerList = header.rstrip("\n").split("\t")
    chromIndex = headerList.index("chrom")
    startIndex = headerList.index("start")
    vcfIDIndex = headerList.index("vcf_id")
    refIndex = headerList.index("ref")
    altIndex = headerList.index("alt")
    geneIndex = headerList.index("gene")
    mscFile.write("chromosome\tposition\tID\treference_allele\talternative_allele\tgene\n")
    for line in chFile:
        lineList = line.rstrip("\n").split("\t")
        chrom = lineList[chromIndex]
        pos = int(lineList[startIndex]) + 1
        vcfID = lineList[vcfIDIndex]
        if vcfID == "None":
            vcfID = "."
        ref = lineList[refIndex]
        alt = lineList[altIndex]
        gene = lineList[geneIndex]
        mscFile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, pos, vcfID, ref, alt, gene))
