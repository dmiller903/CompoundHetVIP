import gzip
import re
import os
import time
import concurrent.futures
import argparse

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='gVCF files are different from VCF files in that they contain information \
for every position, including non-variant positions. Therefore, these files are large and would take a long time to \
process if the non-variant positions were included throughout the whole compound heterozygous pipeline. "keep_variant_sites.py" takes \
each proband file and removes all non-variant sites. If parental files are included, the script then filters each parent \
file for sites that only occur in the proband of that family.')

parser.add_argument('proband_file', help='Proband File')
parser.add_argument('output_path', help='Path to where output files should go')
parser.add_argument('--parent_1_file', help='Maternal or Paternal File of Proband')
parser.add_argument('--parent_2_file', help='Maternal or Paternal File of Proband')
parser.add_argument('--output_suffix', help='Suffix for each output file (do not include .gz at end as this will be included \
when the file is bgzipped', default='_variant_sites.vcf')

args = parser.parse_args()

#Create variables of each argument from argparse
probandFile = args.proband_file
parent1File = args.parent_1_file
parent2File = args.parent_2_file
outputPath = args.output_path
if outputPath.endswith("/"):
    outputPath = outputPath[0:-1]
outputSuffix = args.output_suffix

#Filter each proband file, remove  variants-only sites, create a dictionary of variant-only sites
fileName = re.findall(r'/?([\w\-_]+)\.g\.vcf\.gz', probandFile)[0]
outputName = "{}/{}{}".format(outputPath, fileName, outputSuffix)
positionDict = {}
with gzip.open(probandFile, 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
    for line in gVCF:
        if line.startswith('#'):
            parsed.write(line.encode())
        elif "END" not in line:
            parsed.write(line.encode())
            line = line.split("\t")
            chrom = line[0]
            pos = line[1]
            if chrom not in positionDict:
                positionDict[chrom] = {pos}
            else:
                positionDict[chrom].add(pos)
#bgzip file
os.system("zcat {} | /root/miniconda2/bin/bgzip > {}.gz".format(outputName, outputName))
os.system(f"rm {outputName}")

#Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('Non-variant sites have been removed from probands. Time elapsed: {} minutes ({} hours)'.format(timeElapsedMinutes, timeElapsedHours))

#Filter each parent file for sites that occur in proband of that family
def filterParents(file):
    fileName = re.findall(r'/?([\w\-_]+)\.g\.vcf\.gz', file)[0]
    outputName = "{}/{}{}".format(outputPath, fileName, outputSuffix)
    with gzip.open(file, 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
        for line in gVCF:
            if line.startswith("#"):
                parsed.write(line.encode())
            else:
                lineList = line.split("\t")
                chrom = lineList[0]
                pos = lineList[1]
                if pos in positionDict[chrom]:
                    parsed.write(line.encode())
                else:
                    if "END" in line:
                        for i in range(int(pos), int(lineList[7].lstrip("END=")) + 1):
                            if str(i) in positionDict[chrom]:
                                parsed.write(line.encode())
    #bgzip file
    os.system("zcat {} | /root/miniconda2/bin/bgzip > {}.gz".format(outputName, outputName))
    os.system(f"rm {outputName}")

with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    executor.map(filterParents, [parent1File, parent2File])

#Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('Sites not corresponding to proband have been removed for each parent. Time elapsed: {} minutes ({} hours)'.format(timeElapsedMinutes, timeElapsedHours))