import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Annotate VCF using snpEff")

parser.add_argument('input_vcf', help='Input file')
parser.add_argument('output_file', help='Name of output file')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file

# Annotate the vt trimmed file
os.system(f"java -Xmx40g -jar /snpEff/snpEff.jar GRCh37.75 -v \
{inputFile} > \
{outputFile}")

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Annotation complete. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))