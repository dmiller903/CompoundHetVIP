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
#Separate combinedTrio by chromosome
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

with concurrent.futures.ProcessPoolExecutor(max_workers=42) as executor:
    tmpFileList = executor.map(separateByChr, fileSet)
    for tmpList in tmpFileList:
        for tmpSet in tmpList:
            plinkFileSet.add(tmpSet)

fileDict = dict()
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
                    trioFileName = "{}/{}/{}_trio/{}_trio_{}.vcf".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId, chromosome)
                    #fileDict.add(individualFileName)
                    fileDict[sampleFamilyId].add(trioFileName)

#Create harmonized bed, bim files for each chromosome
for trio in fileDict:
    trioChr = fileDict[trio]
    def createPlink(file):
        filePath, famFolder, sampleFolder, sampleFile, chrNumber = re.findall(r"([\w\-/_]+)\/([\w\-]+)\/([\w\-]+)\/([\w\-]+)_(chr[\w]+)\.?.*\.?.*\.vcf", file)[0]
        outputName = "{}/{}/{}/{}_{}".format(filePath, famFolder, sampleFolder, sampleFile, chrNumber)
        
        os.system("/plink2 --vcf {} --fam {}/{}/{}_trio.fam --make-bed --out {}".format(file, filePath, famFolder, famFolder, outputName))

        if chrNumber[-1].isnumeric():
            os.system("java -jar -Xmx40g /GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar --inputType PLINK_BED --input {} \
            --update-id \
            -ura \
            --outputType PLINK_BED --output {}_harmonized \
            --refType VCF --ref /ALL.{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes".format(outputName, outputName, chrNumber))
        elif chrNumber[-1] == "X":
            os.system("java -jar -Xmx40g /GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar --inputType PLINK_BED --input {} \
            --update-id \
            -ura \
            --outputType PLINK_BED --output {}_harmonized \
            --refType VCF --ref /ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes".format(outputName, outputName, chrNumber))
        else:
            os.system("java -jar -Xmx40g /GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar --inputType PLINK_BED --input {} \
            --update-id \
            -ura \
            --outputType PLINK_BED --output {}_harmonized \
            --refType VCF --ref /ALL.chrY.phase3_integrated_v2a.20130502.genotypes".format(outputName, outputName, chrNumber))

    with concurrent.futures.ProcessPoolExecutor(max_workers=23) as executor:
        executor.map(createPlink, trioChr)

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Trios have been separated by chromosome. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))
