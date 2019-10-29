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

#Create a list of file(s) to be liftedover
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
        sampleIdIndex = header.index("sample_id")
        for sample in sampleFile:
            sampleData = sample.rstrip("\n").split("\t")
            sampleFamilyId = sampleData[familyIdIndex]
            sampleId = sampleData[sampleIdIndex]
            individualFileName = "{}/{}/{}/{}_parsed.vcf.gz".format(pathToFiles, sampleFamilyId, sampleId, sampleId)
            trioFileName = "{}/{}/{}_trio/{}_trio.vcf.gz".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId)
            fileSet.add(individualFileName)
            fileSet.add(trioFileName)

#Liftover file(s)
def liftoverFiles(file):
    if "FM_" in file and "parsed" in file:
        fileFolder, fileName = re.findall("([\w\-\/_]+)\/([\w\-_]+)_parsed\.?.*\.?.*\.gz", file)[0]
        os.system("gatk IndexFeatureFile -F {}".format(file))
        os.system('gatk --java-options "-Xmx4g" GenotypeGVCFs -R /references/Homo_sapiens_assembly38.fasta -V {} -O /tmp/{}_genotyped.vcf.gz'.format(file, fileName))
        os.system("java -jar /root/miniconda2/share/picard-2.21.1-0/picard.jar LiftoverVcf I=/tmp/{}_genotyped.vcf.gz O={}/{}_liftover.vcf.gz \
        CHAIN=/references/hg38ToHg19.over.chain R=/references/human_g1k_v37_modified.fasta REJECT={}/{}_rejected_variants.vcf".format(fileName, fileFolder, fileName, fileFolder, fileName))
    elif "FM_" in file and "parsed" not in file:
        fileFolder, fileName = re.findall("([\w\-\/_]+)\/([\w\-_]+)\.?.*\.?.*\.gz", file)[0]
        os.system("gatk IndexFeatureFile -F {}".format(file))
        os.system('gatk --java-options "-Xmx4g" GenotypeGVCFs -R /references/Homo_sapiens_assembly38.fasta -V {} -O /tmp/{}_genotyped.vcf.gz'.format(file, fileName))
        os.system("java -jar /root/miniconda2/share/picard-2.21.1-0/picard.jar LiftoverVcf I=/tmp/{}_genotyped.vcf.gz O={}/{}_liftover.vcf.gz \
        CHAIN=/references/hg38ToHg19.over.chain R=/references/human_g1k_v37_modified.fasta REJECT={}/{}_rejected_variants.vcf".format(fileName, fileFolder, fileName, fileFolder, fileName))
    else:
        fileName = re.findall("([\w\-\/_\.]+)\.?.*\.?.*\.gz", file)[0]
        os.system("gatk IndexFeatureFile -F {}".format(file))
        os.system('gatk --java-options "-Xmx4g" GenotypeGVCFs -R /references/Homo_sapiens_assembly38.fasta -V {} -O /tmp/{}_genotyped.vcf'.format(file, fileName))
        os.system("java -jar /root/miniconda2/share/picard-2.21.1-0/picard.jar LiftoverVcf I=/tmp/{}_genotyped.vcf O={}_liftover.vcf.gz \
        CHAIN=/references/hg38ToHg19.over.chain R=/references/human_g1k_v37_modified.fasta REJECT={}_rejected_variants.vcf".format(fileName, fileName, fileName))


with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    executor.map(liftoverFiles, fileSet)

#Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Liftover Complete. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))