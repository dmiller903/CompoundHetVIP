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
pathToFiles = argv[2]
numCores = int(argv[3])
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

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
        sampleIdIndex = headerList.index("sample_id")
        for sample in sampleFile:
            sampleData = sample.rstrip("\n").split("\t")
            fileName = sampleData[fileNameIndex]
            sampleFamilyId = sampleData[familyIdIndex]
            sampleId = sampleData[sampleIdIndex]
            individualFileName = "{}/{}/{}/{}_liftover_parsed.vcf.gz".format(pathToFiles, sampleFamilyId, sampleId, sampleId)
            trioFileName = "{}/{}/{}_trio/{}_trio_liftover_parsed.vcf.gz".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId)
            fileSet.add(individualFileName)
            fileSet.add(trioFileName)

plinkFileSet = set()
#Separate combined trio files and individual participant files by chromosome
def separateByChr(file):
    with gzip.open(file, "rt") as vcf:
        fileName = re.findall(r"([\w\-/_]+)_liftover_parsed\.?.*\.?.*\.gz", file)[0]
        outputName = "{}_".format(fileName)
        header = ""
        tmpFileSet = set()
        chromosomeSet = set()
        chromosomeNumber = ""
        
        for line in vcf:
            if line.startswith("#"):
                header = header + line
            elif not line.startswith("#") and line.split("\t")[0] not in chromosomeSet:
                chromosomeNumber = line.split("\t")[0]
                os.system("rm {}{}.*".format(outputName, chromosomeNumber))
                os.system("rm {}{}*harmonized.*".format(outputName, chromosomeNumber))
                os.system("rm {}{}*harmonized*.*".format(outputName, chromosomeNumber))
                with open("{}{}.vcf".format(outputName, chromosomeNumber), "w") as chromosome:
                    chromosome.write(header)
                    chromosome.write(line)
                    chromosomeSet.add(chromosomeNumber)
                    if "trio" in outputName:
                        tmpFileSet.add("{}{}.vcf".format(outputName, chromosomeNumber))
            else:
                with open("{}{}.vcf".format(outputName, chromosomeNumber), "a") as chromosome:
                    chromosome.write(line)
        
        return(tmpFileSet)

with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    tmpFileList = executor.map(separateByChr, fileSet)
    for tmpList in tmpFileList:
        for tmpSet in tmpList:
            plinkFileSet.add(tmpSet)

plinkFileList = list(plinkFileSet)
plinkFileList.sort()

#Create bed, bim files for each chromosome of each trio

def createPlink(file):
    filePath, famFolder, sampleFolder, sampleFile, chrNumber = re.findall(r"([\w\-/_]+)\/([\w\-]+)\/([\w\-]+)\/([\w\-]+)_(chr[\w]+)\.?.*\.?.*\.vcf", file)[0]
    outputName = "{}/{}/{}/{}_{}".format(filePath, famFolder, sampleFolder, sampleFile, chrNumber)

    os.system("/plink2 --vcf {} --fam {}/{}/{}_trio.fam --make-bed --out {}".format(file, filePath, famFolder, famFolder, outputName))

with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    executor.map(createPlink, plinkFileList)


timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Trios have been separated by chromosome. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))