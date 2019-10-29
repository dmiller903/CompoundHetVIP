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
            individualFileName = "{}_test/{}_liftover_parsed.vcf.gz".format(sampleFamilyId, shortName)
            trioFileName = "{}_test/{}_liftover_parsed.vcf.gz".format(sampleFamilyId, sampleFamilyId)
            fileSet.add(individualFileName)
            fileSet.add(trioFileName)

plinkFileSet = set()
#Separate combinedTrio by chromosome
def separateByChr(file):
    with gzip.open(file, "rt") as vcf:
        fileFolder, fileName = re.findall("([\w\-]+)/([\w\-]+)_liftover_parsed\.?.*\.?.*\.gz", file)[0]
        outputName = "{}/{}_".format(fileFolder, fileName)
        header = ""
        tmpFileSet = set()
        for line in vcf:
            if line.startswith("#"):
                header = header + line
            else:
                chromosomeNumber = line.split("\t")[0]
                if not os.path.exists("{}{}.vcf".format(outputName, chromosomeNumber)):
                    with open("{}{}.vcf".format(outputName, chromosomeNumber), "w") as chromosome:
                        chromosome.write(header)
                        chromosome.write(line)
                        if "FM_" in fileName and "{}{}.vcf".format(outputName, chromosomeNumber) not in plinkFileSet:
                            tmpFileSet.add("{}{}.vcf".format(outputName, chromosomeNumber))
                else:
                    with open("{}{}.vcf".format(outputName, chromosomeNumber), "a") as chromosome:
                        chromosome.write(line)
        return(tmpFileSet)

with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    tmpFileList = executor.map(separateByChr, fileSet)
    for tmpList in tmpFileList:
        for tmpSet in tmpList:
            plinkFileSet.add(tmpSet)

print(plinkFileSet)
#Create bed, bim files for each chromosome
def createPlink(file):
    fileFolder, fileName, chrNumber = re.findall("([\w\-]+)\/([\w\-]+)_(chr[\w]+)\.?.*\.?.*\.vcf", file)[0]
    outputName = "{}/{}_{}".format(fileFolder, fileName, chrNumber)
    os.system("/plink2 --vcf {} --fam {}/{}.fam --make-bed --out {}".format(file, fileFolder, fileName, outputName))

with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    executor.map(createPlink, plinkFileSet)

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Trios have been separated by chromosome. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))