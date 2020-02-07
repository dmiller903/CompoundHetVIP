import gzip
import re
import os
import time
from sys import argv
import concurrent.futures
import subprocess

startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

#Combine Manifest and Biospecimen
inputFile = argv[1]
pathToFiles = argv[2]
numCores = int(argv[3])
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

#Create a list of proband files that need to have non-variant sites removed. Create a list of parent files that need sites removed
probandList = []
probandDict = {}
parentList = []
parentDict = {}
with open(inputFile) as tsvFile:
    header = tsvFile.readline()
    header = header.rstrip().split("\t")
    fileNameIndex = header.index("file_name")
    familyIdIndex = header.index("family_id")
    probandIndex = header.index("proband")
    sampleIdIndex = header.index("sample_id")
    for sample in tsvFile:
        sample = sample.rstrip().split("\t")
        if os.path.exists(f"{pathToFiles}/{sample[fileNameIndex]}"):
            if sample[probandIndex] == "Yes":
                probandList.append(sample[fileNameIndex])
                probandDict[sample[fileNameIndex]] = [sample[familyIdIndex], sample[sampleIdIndex]]
            else:
                parentList.append(sample[fileNameIndex])
                parentDict[sample[fileNameIndex]] = [sample[familyIdIndex], sample[sampleIdIndex]]

#Alter parent and proband lists to include only the available downloaded files
newProbandList = []
newParentList = []
for proband in probandList:
    familyName = probandDict[proband][0]
    parentFilesAvailable = 0
    for parentFile in parentDict:
        if parentDict[parentFile][0] == familyName:
            parentFilesAvailable += 1
    if parentFilesAvailable == 2:
        newProbandList.append(proband)
        for parentFile in parentDict:
            if parentDict[parentFile][0] == familyName:
                newParentList.append(parentFile)

#Filter each proband file, remove  variants-only sites, create a dictionary of variant-only sites
positionDict = {}
def filterVariantOnly(file):
    fileName = re.findall(r'(.+)\.g\.vcf\.gz', file)[0]
    familyName = probandDict[file][0]
    sampleName = probandDict[file][1]
    os.system("mkdir {}/{}".format(pathToFiles, familyName))
    os.system("mkdir {}/{}/{}".format(pathToFiles, familyName, sampleName))
    familyDict = {familyName: {}}
    outputName = "{}/{}/{}/{}_parsed.vcf".format(pathToFiles, familyName, sampleName, sampleName)
    with gzip.open("{}/{}".format(pathToFiles, file), 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
        for line in gVCF:
            if line.startswith('#'):
                parsed.write(line.encode())
            elif "END" not in line:
                parsed.write(line.encode())
                line = line.split("\t")
                chrom = line[0]
                pos = line[1]
                if chrom not in familyDict[familyName]:
                    familyDict[familyName][chrom] = {pos}
                else:
                    familyDict[familyName][chrom].add(pos)
    os.system(f"rm {pathToFiles}/{file}")
    return(familyDict)

for i in range(0, len(newProbandList), numCores):
    probandListSlice = newProbandList[i:(i+numCores)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
        executor.map(filterVariantOnly, probandListSlice)
        familyDict = executor.map(filterVariantOnly, probandListSlice)
        for dict in familyDict:
            positionDict.update(dict)

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('Non-variant sites have been removed from probands. Time elapsed: {} minutes ({} hours)'.format(timeElapsedMinutes, timeElapsedHours))

#Filter each parent file for sites that occur in proband of that family
def filterParents(file):
    fileName = re.findall(r'(.+)\.g\.vcf\.gz', file)[0]
    familyName = parentDict[file][0]
    sampleName = parentDict[file][1]
    os.system("mkdir {}/{}/{}".format(pathToFiles, familyName, sampleName))
    outputName = "{}/{}/{}/{}_parsed.vcf".format(pathToFiles, familyName, sampleName, sampleName)
    with gzip.open("{}/{}".format(pathToFiles, file), 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
        for line in gVCF:
            if line.startswith("#"):
                parsed.write(line.encode())
            else:
                lineList = line.split("\t")
                chrom = lineList[0]
                pos = lineList[1]
                if pos in positionDict[familyName][chrom]:
                    parsed.write(line.encode())
                else:
                    if "END" in line:
                        for i in range(int(pos), int(lineList[7].lstrip("END=")) + 1):
                            if str(i) in positionDict[familyName][chrom]:
                                parsed.write(line.encode())
    os.system(f"rm {pathToFiles}/{file}")

for i in range(0, len(newParentList), numCores):
    parentListSlice = newParentList[i:(i+numCores)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
        executor.map(filterParents, parentListSlice)

#bgzip all parsed files
combinedList = newProbandList + newParentList
probandDict.update(parentDict)

def bgzipFiles(file):
    fileName = re.findall(r'(.+)\.g\.vcf\.gz', file)[0]
    familyName = probandDict[file][0]
    sampleName = probandDict[file][1]

    parsedName = "{}/{}/{}/{}_parsed.vcf".format(pathToFiles, familyName, sampleName, sampleName)
    outputName = "{}.gz".format(parsedName)
    os.system("zcat {} | /root/miniconda2/bin/bgzip > {}".format(parsedName, outputName))
    os.system("rm {}".format(parsedName))

for i in range(0, len(combinedList), numCores):
    combinedListSlice = combinedList[i:(i+numCores)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
        executor.map(bgzipFiles, combinedListSlice)

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('Sites not corresponding to proband have been removed for each parent. Time elapsed: {} minutes ({} hours)'.format(timeElapsedMinutes, timeElapsedHours))