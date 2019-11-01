import glob
import re
import os
import concurrent.futures
import subprocess
from sys import argv

#Input file or list of files
inputFile = argv[1]
pathToFiles = argv[2]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]
#Create a list of file(s) that need to have unplaced and multiallelic sites removed
fileDict = dict()
if inputFile.endswith(".gz"):
    fileDict.add(inputFile)

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
        chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
        for sample in sampleFile:
            sampleData = sample.rstrip("\n").split("\t")
            fileName = sampleData[fileNameIndex]
            sampleFamilyId = sampleData[familyIdIndex]
            sampleId = sampleData[sampleIdIndex]
            if sampleFamilyId not in fileDict:
                fileDict[sampleFamilyId] = list()
                for chromosome in chromosomes:
                    #individualFileName = "{}/{}/{}/{}_{}".format(pathToFiles, sampleFamilyId, sampleId, sampleId, chromosome)
                    trioFileName = "{}/{}/{}_trio/{}_trio_{}_phased.vcf".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId, chromosome)
                    #fileDict.add(individualFileName)
                    fileDict[sampleFamilyId].append(trioFileName)


def concatMerge(trio):
    files = fileDict[trio]
    for index, file in enumerate(files):
        os.system("bgzip {} && tabix -p vcf {}.gz".format(file, file))
        files[index] = "{}.gz".format(file)
    fileName = re.findall(r"([\w\-\/_]+)_chr[A-Z0-9][A-Z0-9]?_phased\.vcf", files[0])[0]
    outputName = "{}_phased_combined.vcf".format(fileName)
    files = " ".join(files)
    os.system("bcftools concat {} -o {}".format(files, outputName))
    
with concurrent.futures.ProcessPoolExecutor(max_workers=23) as executor:
    executor.map(concatMerge, fileDict)