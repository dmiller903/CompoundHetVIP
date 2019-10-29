import gzip
import re
import os
import time
import argparse
import concurrent.futures

startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

#Input file
parser = argparse.ArgumentParser(
    description="Combine trio data into one VCF file"
)
parser.add_argument(
    "inputFile",
    help='A TSV with columns "file_name" and "family_id" listing the files to be combined'
)
args = parser.parse_args()

#Create a dictionary of files that need to be combined into one vcf file
fileDict = {}
with open(args.inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        shortName = re.findall(r"([\w\-/]+)\.?.*\.?.*\.gz", fileName)[0]
        actualFileName = "{}_test/{}_parsed.vcf.gz".format(sampleFamilyId, shortName)
        if sampleFamilyId not in fileDict:
            fileDict[sampleFamilyId] = [actualFileName]
        else:
            fileDict[sampleFamilyId].append(actualFileName)

probandDict = {}
parentDict = {}
with open(args.inputFile) as sampleFile:
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
        shortName = re.findall(r"([\w\-/]+)\.?.*\.?.*\.gz", fileName)[0]
        actualFileName = "{}_test/{}_parsed.vcf.gz".format(sampleFamilyId, shortName)
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
    with open(args.inputFile) as sampleFile:
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
    with open("{}_test/{}.fam".format(familyId, familyId), "w") as outputFile:
        for key, value in sorted(sampleDict.items()):
            outputFile.write(value)

with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    executor.map(createFamFiles, probandDict)

filesToGenotype = []
# Use GATK to combine all trios into one vcf
def combineTrios(trio):
    files = fileDict[trio]
    fileString = ""
    outputName = "{}_test/{}.vcf.gz".format(trio, trio)
    for file in files:
        fileString += "-V {} ".format(file)
        os.system("gatk IndexFeatureFile -F {}".format(file))
    os.system("gatk CombineGVCFs -R /references/Homo_sapiens_assembly38.fasta {} -O {}".format(fileString, outputName))
    return(outputName)
with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    outputName = executor.map(combineTrios, fileDict)
    for file in outputName:
        filesToGenotype.append(file)

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Trios have been combined. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))