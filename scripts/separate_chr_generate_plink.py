# Import necessary modules
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
If no fam file is included, PLINK will output a generic fam file. \
see https://www.cog-genomics.org/plink/2.0/formats#fam for formatting guidelines')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file
famFile = args.fam_file

#Separate combined trio files and individual participant files by chromosome
plinkFileSet = set()
with gzip.open(inputFile, "rt") as vcf:
    outputName = f"{outputFile}_"
    chromosomeSet = set()
    chromosomeNumber = ""
    header = ""
    for line in vcf:
        if line.startswith("#"):
            header = header + line
        elif not line.startswith("#") and line.split("\t")[0] not in chromosomeSet:
            chromosomeNumber = line.split("\t")[0]
            if (chromosomeNumber[3:].isnumeric() and int(chromosomeNumber[3:]) in range(0,23)) or (chromosomeNumber[3:].isalpha() and str(chromosomeNumber[3:]) in ["X", "Y"]):
                with gzip.open(f"{outputName}{chromosomeNumber}.vcf.gz", "wb") as chromosome:
                    chromosome.write(header.encode())
                    chromosome.write(line.encode())
                    chromosomeSet.add(chromosomeNumber)
                    plinkFileSet.add(f"{outputName}{chromosomeNumber}.vcf.gz")
            else:
                break
        else:
            with gzip.open(f"{outputName}{chromosomeNumber}.vcf.gz", "ab") as chromosome:
                chromosome.write(line.encode())

#Create bed, bim files for each chromosome of each trio
plinkFileList = list(plinkFileSet)
plinkFileList.sort()
for file in plinkFileList:
    outputName = re.findall(r"([\w/_\-]+chr[0-9XY][0-9XY]?)", file)[0]
    if famFile != None:
        os.system(f"/plink2 --vcf {file} --fam {famFile} --make-bed --out {outputName}")
    else:
        os.system(f"/plink2 --vcf {file} --make-bed --out {outputName}")

# Output Time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')