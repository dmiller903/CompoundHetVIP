import glob
import re
import os
import concurrent.futures
import subprocess
from sys import argv
import gzip

#Input file or list of files
inputFile = argv[1]
pathToFiles = argv[2]
diseaseName = argv[3]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

#Create a list of file(s) that need to have unplaced and multiallelic sites removed
fileDict = dict()
concatFiles = list()
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
                # concatFileName is what each trio will be named after each trio is combined into one trio
                concatFileName = "{}/{}/{}_trio/{}_trio_phased_combined.vcf.gz".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId)
                concatFiles.append(concatFileName)
                for chromosome in chromosomes:
                    #individualFileName = "{}/{}/{}/{}_{}".format(pathToFiles, sampleFamilyId, sampleId, sampleId, chromosome)
                    trioFileName = "{}/{}/{}_trio/{}_trio_{}_phased_reverted.vcf".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId, chromosome)
                    #fileDict.add(individualFileName)
                    fileDict[sampleFamilyId].append(trioFileName)
print(fileDict)
def concatMerge(trio):
    files = fileDict[trio]

    for index, file in enumerate(files):
        #os.system("gzip -d {}.gz".format(file))
        os.system("bgzip -f {} && tabix -fp vcf {}.gz".format(file, file))
        files[index] = "{}.gz".format(file)


    fileName = re.findall(r"([\w\-\/_]+\/[\w\-_]+)_chr[A-Z0-9][A-Z0-9]?_phased_reverted\.vcf", files[0])[0]
    outputName = "{}_phased_combined.vcf".format(fileName)
    files = " ".join(files)
    os.system("bcftools concat {} -o {}".format(files, outputName))
    os.system("bgzip -f {} && tabix -fp vcf {}.gz".format(outputName, outputName))

with concurrent.futures.ProcessPoolExecutor(max_workers=35) as executor:
    executor.map(concatMerge, fileDict)

# Merge all phased, concatenated, trio files into one    
concatFilesString = " ".join(concatFiles)
outputName = "{}/{}_phased_samples.vcf".format(pathToFiles, diseaseName)
os.system("bcftools merge -m both {} -o {}".format(concatFilesString, outputName))
os.system("bgzip -f {} && tabix -fp vcf {}.gz".format(outputName, outputName))

# Create a merged family file
# create a proband dictionary where the key is the sampleId and the value is the familyId
# also create a parent dictionary where the key is familyId and the value is a dictionary that has a key of the sampleId and value of gender
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

# Create a dictionary where each sample has the rest of the family information needed for the family file
sampleDict = dict()
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
        paternal = ""
        maternal = ""
        if probandStatus == "Yes":
            familyDict = parentDict[sampleFamilyId]
            for key, value in familyDict.items():
                if value == "1":
                    paternal = key
                else:
                    maternal = key
            sampleDict[sampleId] = "{}\t{}\t{}\t{}\t{}\t2\n".format(sampleFamilyId, sampleId, paternal, maternal, gender)
        else:
            sampleDict[sampleId] = "{}\t{}\t0\t0\t{}\t1\n".format(sampleFamilyId, sampleId, gender)
            
# create a sample list in the order of the vcf file
with gzip.open("{}/{}_phased_samples.vcf.gz".format(pathToFiles, diseaseName), "rt") as vcfFile:
    for line in vcfFile:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            sampleList = line.rstrip().split("\t")[9:]
        else:
            break

# use the sample order in the list to output each sample in order as found in the vcf file
with open("{}/{}.fam".format(pathToFiles, diseaseName), "w") as outputFile:
    for sample in sampleList:
        outputFile.write(sampleDict[sample])