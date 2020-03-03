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
parser = argparse.ArgumentParser(description="Positions that are multiallelic or are duplicates are removed because \
programs such as PLINK and SHAPEIT2 can not handle these types of sites. Also, sites that contain any missing genotype \
information (i.e. './.') are removed to improve phasing accuracy, but this option can be set to 'n'.")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_vcf', help='Output VCF file')
parser.add_argument('--remove_unknown_genotypes', help="Remove postions with unknown genotypes. This is useful if \
you are working with trios and you are phasing based on family relationships.", default="y")

args = parser.parse_args()

# Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_vcf
removeUnknownGenotypes = args.remove_unknown_genotypes
tempFile = "/tmp/temp.vcf"
fileWithoutSuffix = re.findall(r'([\w\-_/]+)\.', outputFile)[0]
duplicateFile = f"{fileWithoutSuffix}_removed_duplicates.vcf"

# Remove multiallelic sites and positions with unknown genotypes
if removeUnknownGenotypes == "y":
    with gzip.open(inputFile, "rt") as inFile, open(tempFile, "wt") as outFile:
        for line in inFile:
            if line.startswith("#"):
                outFile.write(line) 
            else:
                splitLine = line.split("\t")
                # Only keep positions that are not multiallelic and lines with unknown genotypes
                if "," not in splitLine[4] and "./." not in line:
                    outFile.write(line)
    os.system(f"bgzip -f {tempFile}")
    tempFile = "/tmp/temp.vcf.gz"

# Remove multiallelic sites
else:
    with gzip.open(inputFile, "rt") as inFile, open(tempFile, "wt") as outFile:
        for line in inFile:
            if line.startswith("#"):
                outFile.write(line) 
            else:
                splitLine = line.split("\t")
                # Only keep positions that are not multiallelic
                if "," not in splitLine[4]:
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