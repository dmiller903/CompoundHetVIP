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
parser = argparse.ArgumentParser(description="The parameters for phasing are set so that Eagle uses a haplotype reference panel")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_file', help='Name of output VCF file (without suffix)')
parser.add_argument('chromosome_number', help='Eagle phases chromosome by chromsome, so it needs to be known which \
chromosome to phase')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file
chromosome = args.chromosome_number
tempFile = "/tmp/" + re.findall(r'/?([\w\-_]+)$', outputFile)[-1] + ".vcf"


# VCF files must first have chr# changed to # only
with open(inputFile) as inputFile, open(tempFile, 'w') as output:
    for line in inputFile:
        line = line.replace(f"chr{chromosome}", f"{chromosome}")
        output.write(line)
# Updated VCF needs to be bgzipped and tabixed
os.system("/root/miniconda2/bin/bgzip {}".format(tempFile))
os.system("/root/miniconda2/bin/tabix {}.gz".format(tempFile))
# Phase with Eagle
os.system(f"/Eagle_v2.4.1/eagle --vcfTarget {tempFile}.gz \
--outPrefix {outputFile} \
--geneticMapFile /references/1000GP_Phase3/genetic_map_chr{chromosome}_combined_b37_eagle.txt \
--vcfRef /references/1000GP_Phase3/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")