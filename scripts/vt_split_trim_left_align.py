import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Trim and normalize VCF file")

parser.add_argument('input_vcf', help='Input file')
parser.add_argument('output_file', help='Name of output file')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file

# Use VT to split, trim and left align the phased samples.
os.system(f"/root/miniconda2/bin/vt decompose -s {inputFile} \
| /root/miniconda2/bin/vt normalize -n -r /references/human_g1k_v37.fasta - > \
{outputFile}")

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Trim and normalization complete. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))