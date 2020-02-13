import os
import time
import argparse
import re
import gzip

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Phasing programs require that chromosomes be phased separately. Some \
phasing programs, such as SHAPEIT2, require PLINK files in order to phase. Therefore, this script separates a VCF into \
chromosome VCF files. This step also generates the necessary PLINK files needed for phasing.")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_vcf', help='Path and name of output VCF file')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_vcf
tempFile = "/tmp/temp.vcf"
fileWithoutSuffix = re.findall(r'([\w\-_/]+)\.', outputFile)[0]
duplicateFile = f"{fileWithoutSuffix}_removed_duplicates.vcf"

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