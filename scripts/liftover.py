import os
import time
import argparse
import re

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Some phasing programs and other programs such as GEMINI, require that VCF \
files be aligned to GRCh37. Therefore, this step takes an input VCF that is aligned to GRCh38 converts it to GRCh37 \
positions file using GATKs "CombineGVCFs" tool')

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_vcf', help='Path and name of output VCF file. Do not include ".gz" at the end of the file \
name as this will be included during liftover')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_vcf
fileWithoutSuffix = re.findall(r'([\w\-_/]+)\.', outputFile)[0]

#Liftover file(s)
os.system("java -jar /root/miniconda2/share/picard-2.21.1-0/picard.jar LiftoverVcf I={} O={}.gz \
CHAIN=/references/hg38ToHg19.over.chain R=/references/human_g1k_v37_modified.fasta REJECT={}_rejected_variants.vcf".format(inputFile, outputFile, fileWithoutSuffix))

#Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Liftover Complete. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))