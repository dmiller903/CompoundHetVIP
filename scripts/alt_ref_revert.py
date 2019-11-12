import glob
import gzip
import re
from sys import argv

inputFile = argv[1]
pathToFiles = argv[2]

fileDict = dict()
if inputFile.endswith(".vcf"):
    fileDict["1"] = {inputFile}

elif inputFile.endswith(".txt"):
    with open(inputFile) as sampleFile:
        for sample in sampleFile:
            sample = sample.rstrip("\n")
            fileDict.add(sample)

elif inputFile.endswith(".tsv"):
    with open(inputFile) as sampleFile:
            header = sampleFile.readline()
            headerList = header.rstrip().split("\t")
            fileNameIndex = headerList.index("file_name")
            familyIdIndex = headerList.index("family_id")
            sampleIdIndex = headerList.index("sample_id")
            chromosomes = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
    "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"}
            for sample in sampleFile:
                sampleData = sample.rstrip("\n").split("\t")
                fileName = sampleData[fileNameIndex]
                sampleFamilyId = sampleData[familyIdIndex]
                sampleId = sampleData[sampleIdIndex]
                if sampleFamilyId not in fileDict:
                    fileDict[sampleFamilyId] = set()
                    for chromosome in chromosomes:
                        #individualFileName = "{}/{}/{}/{}_{}".format(pathToFiles, sampleFamilyId, sampleId, sampleId, chromosome)
                        trioFileName = "{}/{}/{}_trio/{}_trio_{}_phased.vcf".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId, chromosome)
                        #fileDict.add(individualFileName)
                        fileDict[sampleFamilyId].add(trioFileName)

posDict = dict()
for file in glob.glob("/references/1000GP_Phase3/*legend.gz"):
    with gzip.open(file, 'rt') as legend:
        chrom = re.findall(r"[\w_/]+_chr([0-9]+)\.legend\.gz", file)[0]
        header = legend.readline()
        headerList = header.rstrip().split()
        refIndex = headerList.index("a0")
        altIndex = headerList.index("a1")
        posIndex = headerList.index("position")
        idIndex = headerList.index("id")
        for line in legend:
            lineList = line.rstrip().split(" ")
            pos = lineList[posIndex]
            ref = lineList[refIndex]
            alt = lineList[altIndex]
            siteStr = "{} {} {}".format(pos, ref, alt)
            if chrom not in posDict:
                posDict[chrom] = {siteStr}
            else:
                posDict[chrom].add(siteStr)

print("Dictionary Created\n")

for key, value in fileDict.items():
    for file in value:
        rawCount = 0
        flipCount = 0
        total = 0
        outputName = re.findall(r"([\w\-\/_]+\/[\w\-_]+_chr[A-Z0-9][A-Z0-9]?[_\w]*_phased)\.vcf", file)[0]
        with open(file, 'rt') as sample, open("{}_reverted.vcf".format(outputName), 'w') as output:
            for line in sample:
                if "##" in line:
                    output.write(line)
                elif line.startswith("#CHROM"):
                    header = line.split("\t")
                    chromIndex = header.index("#CHROM")
                    posIndex = header.index("POS")
                    refIndex = header.index("REF")
                    altIndex = header.index("ALT")
                    output.write(line)
                else:
                    lineList = line.split("\t")
                    chrom = lineList[chromIndex]
                    pos = lineList[posIndex]
                    ref = lineList[refIndex]
                    alt = lineList[altIndex]
                    rawStr = "{} {} {}".format(pos, ref, alt)
                    flipStr = "{} {} {}".format(pos, alt, ref)
                    if rawStr in posDict[chrom]:
                        output.write(line)
                        rawCount += 1
                        total += 1
                    elif flipStr in posDict[chrom]:
                        lineList[refIndex] = alt
                        lineList[altIndex] = ref
                        line = "\t".join(lineList)
                        line = line.replace("0|1", "b|a").replace("1|0", "a|b").replace("1|1", "a|a").replace("0|0", "b|b")
                        line = line.replace("b|a", "1|0").replace("a|b", "0|1").replace("a|a", "0|0").replace("b|b", "1|1")
                        output.write(line)
                        flipCount += 1
                        total += 1
                    else:
                        total += 1
            print(file)
            rawPercent = (rawCount / total) * 100
            flipPercent = (flipCount / total) * 100
            totalPercent = ((flipCount + rawCount) / total) * 100
            print("For {}, chr{}, {} ({:.2f}%) of the sites were unchanged".format(key, chrom, rawCount, rawPercent))
            print("For {}, chr{}, {} ({:.2f}%) of the sites were switched to match the reference panel".format(key, chrom, flipCount, flipPercent))
            print("For {}, chr{}, {:.2f}% of the sites are now congruent with the reference panel\n".format(key, chrom, totalPercent))