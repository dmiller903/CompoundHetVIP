#  Import necessary modules
import os
import time
import argparse
import glob

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Creates a combined trio file using GATKs "CombineGVCFs" and GenotypeGVCFs \
tools if gVCF files are used. If VCF files are used, the trio will be combined with bcftools')

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

# Download reference files if needed
if not os.path.exists("/references/Homo_sapiens_assembly38.fasta") and isGvcf == "y":
    os.system("wget --no-check-certificate \
    https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d9f54d2a7bc73000ee99fd6/?zip= -O /tmp/references.zip \
    && unzip /tmp/references.zip -d /tmp/references \
    && rm /tmp/references.zip \
    && gzip -d /tmp/references/*.gz")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}"):
            os.system(f"mv {file} /references/")
    os.system("chmod 777 /references/*")

# Use GATK to combine all trios into one vcf and then genotype the combined trio vcf
files = [probandFile, parent1File, parent2File]
tempName = "/tmp/temp.vcf.gz"
if isGvcf == "y":
    try:
        fileString = ""
        for file in files:
            fileString += f"-V {file} "
            os.system(f"/root/miniconda2/bin/gatk IndexFeatureFile -F {file}")
        os.system(f"/root/miniconda2/bin/gatk CombineGVCFs -R /references/Homo_sapiens_assembly38.fasta {fileString} -O {tempName}")
        os.system(f"gatk IndexFeatureFile -F {tempName}")
        os.system(f"gatk --java-options '-Xmx4g' GenotypeGVCFs -R /references/Homo_sapiens_assembly38.fasta -V {tempName} -O {outputName}")
    except:
        print("Trio not combined, there was an error detected by GATK")

elif isGvcf == "n":
    for file in files:
        os.system(f"tabix -fp vcf {file}")
    fileString = " ".join(files)
    os.system(f"bcftools merge {fileString} -o {outputName.rstrip(".gz")} && bgzip {outputName.rstrip(".gz")} && tabix -fp vcf {outputName}")

# Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')