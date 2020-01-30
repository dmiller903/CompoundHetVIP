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

for trio in fileDict:
    trioChr = fileDict[trio]
    def runEagle(file):
        chromosome = re.findall(r"[\w\-\/_]+chr([0-9][0-9]?)", file)[0]
        outputName = "/tmp/" + re.findall(r"\/([\w\-_]+chr[0-9][0-9]?)", file)[0] + "_update.vcf"
        filePrefix = re.findall(r"([\w\-\/_]+chr[0-9][0-9]?)", file)[0] + "_eagle_phased"
        # VCF files must first have chr# changed to # only
        with open(file) as inputFile, open(outputName, 'w') as outputFile:
            for line in inputFile:
                line = line.replace(f"chr{chromosome}", f"{chromosome}")
                outputFile.write(line)
        # Updated VCF needs to be bgzipped and tabixed
        os.system("/root/miniconda2/bin/bgzip {}".format(outputName))
        os.system("/root/miniconda2/bin/tabix {}.gz".format(outputName))
        # Phase with Eagle
        os.system(f"java -Xmx40g -jar /beagle.25Nov19.28d.jar gt={outputName}.gz \
        out={filePrefix} \
        chrom={chromosome} \
        map=/references/1000GP_Phase3/genetic_map_chr{chromosome}_combined_b37_eagle.txt \
        ref=/references/1000GP_Phase3/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        impute=false")
    with concurrent.futures.ProcessPoolExecutor(max_workers=23) as executor:
        executor.map(runEagle, trioChr)


