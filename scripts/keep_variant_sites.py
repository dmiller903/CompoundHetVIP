# Import necessary modules
import gzip
import re
import os
import time
import concurrent.futures
import argparse

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='If input files are VCF and parent VCFs are unavailable or not being \
used, this script does not need to be used. Whether using gVCF or VCF files, they must be gzipped. When using VCF \
files as input, the optional parameter "--is_gvcf" needs to be changed to "n". This script will filter each parent \
VCF for sites that only occur in the child of that family. When using gVCF files, this scrip behaves a little \
differently. gVCF files are different from VCF files in that they contain information for every nucleotide position, \
including non-variant positions. Therefore, gVCF files are large and would take a long time to process if the \
non-variant positions were included throughout the whole pipeline. Thus, this script takes the sample (patient) gVCF \
and removes all non-variant sites. In addition, if parental gVCFs are included, the script filters each parent file for \
sites that only occur in the affected sample of that family.')

parser.add_argument('sample_file', help='Sample (patient) File. Must be gzipped')
parser.add_argument('output_path', help='Path to where output file(s) should go')
parser.add_argument('--parent_1_file', help='Maternal or Paternal File of Sample. Must be gzipped')
parser.add_argument('--parent_2_file', help='Maternal or Paternal File of Sample. Must be gzipped')
parser.add_argument('--output_suffix', help='Suffix for each output file (output file will be bgzipped)', default='_parsed.vcf.gz')
parser.add_argument('--is_gvcf', help='If a gVCF file is used, all non-variant sites will be filtered out of the sample \
file and a new VCF will be created', default='y')

args = parser.parse_args()

#Create variables of each argument from argparse
sampleFile = args.sample_file
parent1File = args.parent_1_file
parent2File = args.parent_2_file
outputPath = args.output_path.rstrip("/")
isGvcf = args.is_gvcf
outputSuffix = args.output_suffix.rstrip(".gz")

# Set of chromosomes to keep
chrToKeep = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

#Filter each sample file, remove  variants-only sites, create a dictionary of variant-only sites
fileName = re.findall(r'/?([\w\-_]+)\.?g?\.vcf\.gz', sampleFile)[0]
outputName = f"{outputPath}/{fileName}{outputSuffix}"
positionDict = {}
if isGvcf == "y":
    with gzip.open(sampleFile, 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
        for line in gVCF:
            if line.startswith('#'):
                parsed.write(line.encode())
            elif "END=" not in line:
                line_list = line.split("\t")
                chrom = line_list[0]
                pos = line_list[1]
                if chrom not in positionDict and chrom in chrToKeep:
                    positionDict[chrom] = {pos}
                    parsed.write(line.encode())
                elif chrom in positionDict:
                    positionDict[chrom].add(pos)
                    parsed.write(line.encode())
    #bgzip file
    os.system(f"zcat {outputName} | /root/miniconda2/bin/bgzip > {outputName}.gz")
    os.system(f"rm {outputName}")

elif isGvcf == "n":
    with gzip.open(sampleFile, 'rt') as gVCF:
        for line in gVCF:
            if line.startswith('#'):
                continue
            else:
                line = line.split("\t")
                chrom = line[0]
                pos = line[1]
                if chrom not in positionDict and chrom in chrToKeep:
                    positionDict[chrom] = {pos}
                elif chrom in positionDict:
                    positionDict[chrom].add(pos)

#Filter each parent file for sites that occur in sample of that family
def filterParents(file):
    fileName = re.findall(r'/?([\w\-_]+)\.?g?\.vcf\.gz', file)[0]
    outputName = f"{outputPath}/{fileName}{outputSuffix}"
    with gzip.open(file, 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
        for line in gVCF:
            if line.startswith("#"):
                parsed.write(line.encode())
            else:
                lineList = line.split("\t")
                chrom = lineList[0]
                pos = lineList[1]
                if chrom in positionDict and pos in positionDict[chrom]:
                    parsed.write(line.encode())
                else:
                    if "END=" in line and isGvcf == "y":
                        for i in range(int(pos), int(lineList[7].lstrip("END=")) + 1):
                            if chrom in positionDict and str(i) in positionDict[chrom]:
                                parsed.write(line.encode())
        return(outputName)

if parent1File != None and parent2File != None:
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        for outputName in executor.map(filterParents, [parent1File, parent2File]):
            #bgzip file
            os.system(f"zcat {outputName} | /root/miniconda2/bin/bgzip > {outputName}.gz")
            os.system(f"rm {outputName}")

#Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours) {char}')
