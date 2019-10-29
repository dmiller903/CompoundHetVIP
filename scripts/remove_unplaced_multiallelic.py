import gzip
import re
import os
import time
from sys import argv
import concurrent.futures

startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

#Input file or list of files
inputFile = argv[1]

#Create a list of file(s) that need to have unplaced and multiallelic sites removed
fileSet = set()
if inputFile.endswith(".gz"):
    fileSet.add(inputFile)

elif inputFile.endswith(".txt"):
    with open(inputFile) as sampleFile:
        for sample in sampleFile:
            sample = sample.rstrip("\n")
            fileSet.add(sample)

elif inputFile.endswith(".tsv"):
    with open(inputFile) as sampleFile:
        header = sampleFile.readline()
        headerList = header.rstrip().split("\t")
        fileNameIndex = headerList.index("file_name")
        familyIdIndex = headerList.index("family_id")
        for sample in sampleFile:
            sampleData = sample.rstrip("\n").split("\t")
            fileName = sampleData[fileNameIndex]
            sampleFamilyId = sampleData[familyIdIndex]
            shortName = re.findall(r"([\w\-/]+)\.?.*\.?.*\.gz", fileName)[0]
            individualFileName = "{}_test/{}_liftover.vcf.gz".format(sampleFamilyId, shortName)
            trioFileName = "{}_test/{}_liftover.vcf.gz".format(sampleFamilyId, sampleFamilyId)
            fileSet.add(individualFileName)
            fileSet.add(trioFileName)

#Remove Unplaced sites and multiallelic sites
chrToKeep = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

filesToRemoveDuplicates = []

def removeSites(file):
    fileName = re.findall(r"([\w\-/]+)_liftover\.?.*\.?.*\.gz", file)[0]
    outputName = "{}_no_ambiguous_sites.vcf".format(fileName)
    with gzip.open(file, "rt") as inputFile, open(outputName, "wt") as outFile:
        for line in inputFile:
            if line.startswith("#") and "##contig=<ID=" not in line:
                outFile.write(line)   
            elif line.startswith("#") and "##contig=<ID=" in line:
                splitLine = line.split(",")
                chr = splitLine[0].replace("##contig=<ID=", "")
                if chr in chrToKeep:
                    outFile.write(line)
            else:
                splitLine = line.split("\t")
                if splitLine[0] in chrToKeep and "," not in splitLine[4]:
                    outFile.write(line)
    os.system("bgzip {}".format(outputName))
    return("{}.gz".format(outputName))

with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    fileName = executor.map(removeSites, fileSet)
    for file in fileName:
        filesToRemoveDuplicates.append(file)

#Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Unplaced and multiallelic sites have been removed. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))

#Remove all duplicate sites
def removeDuplicates(file):
    fileName = re.findall(r"([\w\-/]+)_no_ambiguous_sites\.?.*\.?.*\.gz", file)[0]
    outputName = "{}_liftover_parsed.vcf".format(fileName)
    duplicateFile = "{}_removedDuplicates.vcf".format(fileName)

    lines = [["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"], ["NA"]]

    dupList = []
    with gzip.open(file, "rt") as inputFile:
        for line in inputFile:
            if not line.startswith("#"):
                line = line.split("\t")
                line = line[0:2] + line[3:5]
                if line not in lines:
                    lines.append(line)
                    lines = lines[-10:]
                else:
                    dupList.append(line)

    with gzip.open(file, "rt") as inputFile, open(outputName, "wt") as outFile, open(duplicateFile, "w") as duplicates:
        for line in inputFile:
            if not line.startswith("#"):
                splitLine = line.split("\t")
                splitLine = splitLine[0:2] + splitLine[3:5]
                if splitLine not in dupList:
                    outFile.write(line)
                else:
                    duplicates.write(line)
            else:
                outFile.write(line)
                duplicates.write(line)
    os.system("bgzip {}".format(outputName))
with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    executor.map(removeDuplicates, filesToRemoveDuplicates)


timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Duplicate sites removed. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))
