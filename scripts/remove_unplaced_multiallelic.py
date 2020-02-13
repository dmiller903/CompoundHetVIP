import os
import time
import argparse
import re
import gzip

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="During liftover, some randomly placed sites are included in the VCF file. \
These randomly placed sites are those that are in GRCh38 but the exact position in GRCh37 isn't known. Therefore, for \
subsequent analysis, these sites are removed. Only sites with known positions, on a known chromosome are kept. In \
addition, positions that are multiallelic or are duplicates are removed because programs such as PLINK and \
SHAPEIT2 can not handle these types of sites. Also, sites that contain any missing genotype information \
(i.e. './.') are removed to improve phasing accuracy")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_vcf', help='Path and name of output VCF file')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_vcf
tempFile = "/tmp/temp.vcf"
fileWithoutSuffix = re.findall(r'([\w\-_/]+)\.', outputFile)[0]
duplicateFile = f"{fileWithoutSuffix}_removed_duplicates.vcf"

#Remove Unplaced sites, multiallelic sites, and sites where the genotype ./. occurs more than once
chrToKeep = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

with gzip.open(inputFile, "rt") as inFile, open(tempFile, "wt") as outFile:
    for line in inFile:
        if line.startswith("#") and "##contig=<ID=" not in line:
            outFile.write(line)   
        elif line.startswith("#") and "##contig=<ID=" in line:
            splitLine = line.split(",")
            chr = splitLine[0].replace("##contig=<ID=", "")
            if chr in chrToKeep:
                outFile.write(line)
        else:
            splitLine = line.split("\t")
            # Only keep positions where genotype information is available for all samples
            if splitLine[0] in chrToKeep and "," not in splitLine[4] and "./." not in line:
                outFile.write(line)

os.system("bgzip -f {}".format(tempFile))
tempFile = "/tmp/temp.vcf.gz"

#Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Unplaced sites, multiallelic sites, and sites where ./. occurs more than once have been removed. Time elapsed: \
{} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))

#Remove all duplicate sites
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

os.system("bgzip -f {}".format(outputFile))


timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Duplicate sites removed. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))
