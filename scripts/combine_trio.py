import os
import time
import argparse

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Creates a combined trio file using GATKs "CombineGVCFs" and GenotypeGVCFs \
tools')

parser.add_argument('proband_vcf', help='Proband VCF File')
parser.add_argument('parent_1_vcf', help='Maternal or Paternal VCF File of Proband')
parser.add_argument('parent_2_vcf', help='Maternal or Paternal VCF File of Proband')
parser.add_argument('output_vcf', help='Path and name of combined vcf output file')
parser.add_argument('--is_gvcf', help="If a gVCF file is used, GATK's combineGVCFs will be used. If VCF files are used \
bcftools merge will be used", default='y')

args = parser.parse_args()

#Create variables of each argument from argparse
probandFile = args.proband_vcf
parent1File = args.parent_1_vcf
parent2File = args.parent_2_vcf
outputName = args.output_vcf
isGvcf = args.is_gvcf

# Use GATK to combine all trios into one vcf and then genotype the combined trio vcf
files = [probandFile, parent1File, parent2File]
tempName = "/tmp/temp.vcf"
if isGvcf == "y":
    fileString = ""
    for file in files:
        fileString += "-V {} ".format(file)
        os.system("/root/miniconda2/bin/gatk IndexFeatureFile -F {}".format(file))
    os.system("/root/miniconda2/bin/gatk CombineGVCFs -R /references/Homo_sapiens_assembly38.fasta {} -O {}".format(fileString, tempName))
    os.system("gatk IndexFeatureFile -F {}".format(tempName))
    os.system('gatk --java-options "-Xmx4g" GenotypeGVCFs -R /references/Homo_sapiens_assembly38.fasta -V {} -O {}'.format(tempName, outputName))

elif isGvcf == "n":
    for file in files:
        os.system("tabix -fp vcf {}".format(file))
    fileString = " ".join(files)
    os.system("bcftools merge {} -o {}".format(fileString, outputName))
    os.system("bgzip -f {} && tabix -fp vcf {}.gz".format(outputName, outputName))

#Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Trio has been combined. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))