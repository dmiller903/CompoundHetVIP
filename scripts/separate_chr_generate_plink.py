import os
import time
import argparse
import re
import gzip
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Phasing programs require that chromosomes be phased separately. Some \
phasing programs, such as SHAPEIT2, require PLINK files in order to phase. Therefore, this script separates a VCF into \
chromosome VCF files. This step also generates the necessary PLINK files needed for phasing.")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_file', help='Path and prefix name of output files (no suffix). e.g. "/Data/file1"')
parser.add_argument('--fam_file', help='If using a trio, a fam file is needed to create appropriate PLINK files. \
If no fam file is included, PLINK with output a generic fam file. \
see https://www.cog-genomics.org/plink/2.0/formats#fam for formatting guidelines')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file
famFile = args.fam_file

#Separate combined trio files and individual participant files by chromosome
plinkFileSet = set()
with gzip.open(inputFile, "rt") as vcf:
    outputName = "{}_".format(outputFile)
    chromosomeSet = set()
    chromosomeNumber = ""
    header = ""
    for line in vcf:
        if line.startswith("#"):
            header = header + line
        elif not line.startswith("#") and line.split("\t")[0] not in chromosomeSet:
            chromosomeNumber = line.split("\t")[0]
            if chromosomeNumber.startswith("c"):
                chromosomeNumber = chromosomeNumber[3:]
                if (chromosomeNumber.isnumeric() and int(chromosomeNumber) in range(0,23)) or (chromosomeNumber.isalpha() and str(chromosomeNumber) in ["X", "Y"]):
                    with open("{}chr{}.vcf".format(outputName, chromosomeNumber), "w") as chromosome:
                        chromosome.write(header)
                        lineList = line.split("\t")
                        updateChromosome = lineList[0]
                        updateChromosome = updateChromosome[3:]
                        lineList[0] = updateChromosome
                        line = " ".join(lineList)
                        chromosome.write(line)
                        chromosomeSet.add("chr" + chromosomeNumber)
                        plinkFileSet.add("{}chr{}.vcf".format(outputName, chromosomeNumber))
                else:
                    break
            else:
                if (chromosomeNumber.isnumeric() and int(chromosomeNumber) in range(0,23)) or (chromosomeNumber.isalpha() and str(chromosomeNumber) in ["X", "Y"]):
                    with open("{}chr{}.vcf".format(outputName, chromosomeNumber), "w") as chromosome:
                        chromosome.write(header)
                        chromosome.write(line)
                        chromosomeSet.add(chromosomeNumber)
                        plinkFileSet.add("{}chr{}.vcf".format(outputName, chromosomeNumber))
                else:
                    break
        else:
            if chromosomeNumber.startswith("c"):
                with open("{}chr{}.vcf".format(outputName, chromosomeNumber), "a") as chromosome:
                    lineList = line.split("\t")
                    updateChromosome = lineList[0]
                    updateChromosome = updateChromosome[3:]
                    lineList[0] = updateChromosome
                    line = " ".join(lineList)
                    chromosome.write(line)
            else:
                with open("{}chr{}.vcf".format(outputName, chromosomeNumber), "a") as chromosome:
                    chromosome.write(line)

#Create bed, bim files for each chromosome of each trio
plinkFileList = list(plinkFileSet)
print(plinkFileList)
plinkFileList.sort()
for file in plinkFileList:
    outputName = re.findall(r"([\w/_\-]+chr[0-9XY][0-9XY]?)", file)[0]
    if famFile != None:
        os.system("/plink2 --vcf {} --fam {} --make-bed --out {}".format(file, famFile, outputName))
    else:
        os.system("/plink2 --vcf {} --make-bed --out {}".format(file, outputName))

timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}File has been separated by chromosome. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))