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

#Create a dictionary of files that need to be combined into one vcf file
fileDict = {}
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        print(sampleData)
        sampleId = sampleData[sampleIdIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        actualFileName = "{}/{}/{}/{}_parsed.vcf.gz".format(pathToFiles, sampleFamilyId, sampleId, sampleId)
        if sampleFamilyId not in fileDict:
            fileDict[sampleFamilyId] = [actualFileName]
        else:
            fileDict[sampleFamilyId].append(actualFileName)

probandDict = {}
parentDict = {}
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    probandIndex = headerList.index("proband")
    genderIndex = headerList.index("sex")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        probandStatus = sampleData[probandIndex]
        gender = sampleData[genderIndex]
        if probandStatus == "Yes":
            probandDict[sampleId] = sampleFamilyId
        else:
            if sampleFamilyId not in parentDict:
                parentDict[sampleFamilyId] = {sampleId: gender}
            else:
                parentDict[sampleFamilyId][sampleId] = gender

# Create fam files
def createFamFiles(proband):
    familyId = probandDict[proband]
    familyDict = parentDict[familyId]
    paternal = ""
    maternal = ""
    outputString = ""
    sampleDict = {}
    for key, value in familyDict.items():
        if value == "1":
            paternal = key
        else:
            maternal = key
    with open(inputFile) as sampleFile:
        header = sampleFile.readline()
        headerList = header.rstrip().split("\t")
        fileNameIndex = headerList.index("file_name")
        familyIdIndex = headerList.index("family_id")
        sampleIdIndex = headerList.index("sample_id")
        probandIndex = headerList.index("proband")
        genderIndex = headerList.index("sex")
        for sample in sampleFile:
            sampleData = sample.rstrip("\n").split("\t")
            fileName = sampleData[fileNameIndex]
            sampleFamilyId = sampleData[familyIdIndex]
            sampleId = sampleData[sampleIdIndex]
            probandStatus = sampleData[probandIndex]
            gender = sampleData[genderIndex]
            if probandStatus == "Yes" and familyId == sampleFamilyId:
                sampleDict[sampleId] = "{}\t{}\t{}\t{}\t{}\t2\n".format(sampleFamilyId, sampleId, paternal, maternal, gender)
            elif probandStatus == "No" and familyId == sampleFamilyId:
                sampleDict[sampleId] = "{}\t{}\t0\t0\t{}\t1\n".format(sampleFamilyId, sampleId, gender)
    with open("{}/{}/{}_trio.fam".format(pathToFiles, familyId, familyId), "w") as outputFile:
        for key, value in sorted(sampleDict.items()):
            outputFile.write(value)

with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    executor.map(createFamFiles, probandDict)

filesToGenotype = []
# Use GATK to combine all trios into one vcf
def combineTrios(trio):
    files = fileDict[trio]
    fileString = ""
    os.system("mkdir {}/{}/{}_trio".format(pathToFiles, trio, trio))
    outputName = "{}/{}/{}_trio/{}_trio.vcf.gz".format(pathToFiles, trio, trio, trio)
    for file in files:
        fileString += "-V {} ".format(file)
        os.system("gatk IndexFeatureFile -F {}".format(file))
    os.system("gatk CombineGVCFs -R /references/Homo_sapiens_assembly38.fasta {} -O {}".format(fileString, outputName))
    return(outputName)
with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    outputName = executor.map(combineTrios, fileDict)
    for file in outputName:
        filesToGenotype.append(file)

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Trios have been combined. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))