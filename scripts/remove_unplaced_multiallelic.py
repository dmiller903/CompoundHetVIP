# import necessary modules
import os
import time
import argparse
import re
import gzip

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Positions that are multiallelic or duplicates are removed because \
programs such as PLINK and SHAPEIT2 can not handle these types of sites. Also, positions that contain missing genotype \
information (i.e. './.') in more then one sample are removed to improve phasing accuracy.")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_vcf', help='Output VCF file')

args = parser.parse_args()

# Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_vcf.rstrip(".gz")
tempFile = "/tmp/temp.vcf"
fileWithoutSuffix = re.findall(r'([\w\-_/]+)\.', outputFile)[0]
duplicateFile = f"{fileWithoutSuffix}_removed_duplicates.vcf"

# Set of chromosomes to keep
chrToKeep = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

# Remove multiallelic sites and keep positions where genotype information is available for patient and at least one parent
with gzip.open(inputFile, "rt") as inFile, open(tempFile, "wt") as outFile:
    for line in inFile:
        if line.startswith("#"):
            outFile.write(line) 
        else:
            splitLine = line.split("\t")
            if splitLine[0] in chrToKeep and "," not in splitLine[4] and line.count("./.") < 2:
                outFile.write(line)
                
os.system(f"bgzip -f {tempFile}")
tempFile = "/tmp/temp.vcf.gz"

# Remove all duplicate sites
posDict = dict()
dupDict = dict()
with gzip.open(tempFile, "rt") as inputFile:
    for line in inputFile:
        if not line.startswith("#"):
            line = line.split("\t")
            chromosome = line[0]
            pos = line[1]
            if chromosome not in posDict:
                posDict[chromosome] = set()
                posDict[chromosome].add(pos)
                dupDict[chromosome] = set()
            elif chromosome in posDict and pos not in posDict[chromosome]:
                posDict[chromosome].add(pos)
            elif chromosome in posDict and pos in posDict[chromosome]:
                dupDict[chromosome].add(pos)

with gzip.open(tempFile, "rt") as inputFile, open(outputFile, "wt") as outFile, open(duplicateFile, "w") as duplicates:
    for line in inputFile:
        if not line.startswith("#"):
            splitLine = line.split("\t")
            chromosome = splitLine[0]
            pos = splitLine[1]
            if pos not in dupDict[chromosome]:
                outFile.write(line)
            else:
                duplicates.write(line)
        else:
            outFile.write(line)
            duplicates.write(line)

os.system(f"bgzip -f {outputFile}")

# Output time it took to complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')