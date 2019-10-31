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
                    trioFileName = "{}/{}/{}_trio/{}_trio_{}".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId, chromosome)
                    #fileDict.add(individualFileName)
                    fileDict[sampleFamilyId].add(trioFileName)

for trio in fileDict:
    trioChr = fileDict[trio]
    def runShapeit(file):
        chromosome = re.findall(r"[\w\-\/_]+(chr[0-9][0-9]?)", file)[0]
        os.system("/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -check -B {} --output-log {} -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3".format(file, file, chromosome, chromosome, chromosome))
        os.system("/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -B {} --output-log {} -O {}_phased -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --exclude-snp {}.snp.strand.exclude --force --no-mcmc".format(file, file, file, chromosome, chromosome, chromosome, file))
        #os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {} -M /references/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_{}_combined_b37.txt \
                #--input-ref /references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3.sample /references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3_{}_impute_macGT1.hap.gz \
                #/references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3_{}_impute_macGT1.legend.gz -O {}_phased --output-graph {}_phased.graph \
                #--output-log {}_phased.log --exclude-snp /fslhome/dmill903/research/idiopathic_scoliosis/gVCF/FM_YTPYPGH1_test/FM_YTPYPGH1_chr21.checks.snp.strand.exclude".format(file, chr, chr, chr, file, file, file))
        os.system("/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -convert --input-haps {}_phased --output-vcf {}_phased.vcf".format(file, file))
    with concurrent.futures.ProcessPoolExecutor(max_workers=23) as executor:
        executor.map(runShapeit, trioChr)
