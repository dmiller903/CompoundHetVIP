# import necessary modules
import os
import time
import argparse
import re
import gzip
import concurrent.futures
import glob

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="The parameters for phasing are set so that Beagle uses a haplotype reference panel")

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_file', help='Name of output VCF file (without suffix)')
parser.add_argument('chromosome_number', help='Beagle phases chromosome by chromsome, so it needs to be known which \
chromosome to phase')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file
chromosome = args.chromosome_number
tempFile = "/tmp/" + re.findall(r'/?([\w\-_]+)$', outputFile)[-1] + ".vcf"

# Download reference files if necessary
if not os.path.exists("/references/1000GP_Phase3/1000GP_Phase3.sample"):
    os.system("wget --no-check-certificate https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz -P /tmp/references \
    && tar zxf /tmp/references/1000GP_Phase3.tgz -C /tmp/references/ \
    && rm /tmp/references/1000GP_Phase3.tgz \
    && wget --no-check-certificate \
    https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5db8af02f3bb87000b85b76e/?zip= -O /tmp/references.zip \
    && unzip /tmp/references.zip -d /tmp/references/1000GP_Phase3 \
    && rm /tmp/references.zip")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}"):
            os.system(f"mv {file} /references/")
    os.system("chmod -R 777 /references/1000GP_Phase3")

# Modify genetic map files to format required for beagle if not already done
if not os.path.exists("/references/1000GP_Phase3/genetic_map_chr1_combined_b37_beagle.txt"):
    def updateFiles(file):
        fileName = re.findall(r"([\w/_]+genetic_map_chr(\w+)_combined_b37)\.txt", file)[0][0]
        chrom = re.findall(r"([\w/_]+genetic_map_chr(\w+)_combined_b37)\.txt", file)[0][1]
        beagleOutput = f"{fileName}_beagle.txt"
        with open(file) as inputFile, open(beagleOutput, 'w') as beagleOut:
            header = inputFile.readline()
            header = "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
            for line in inputFile:
                beagleLineList = line.rstrip().split(" ")
                beagleLine = f"{chrom} {beagleLineList[1]} {beagleLineList[2]} {beagleLineList[0]}\n"
                beagleOut.write(beagleLine)
        os.system(f"chmod 777 {beagleOutput}")
    for file in glob.glob("/references/1000GP_Phase3/genetic_map_chr*_combined_b37.txt"):
        updateFiles(file)

# VCF files must first have chr# changed to # only
with gzip.open(inputFile, "rt") as inputFile, open(tempFile, 'w') as output:
    for line in inputFile:
        line = line.replace(f"chr{chromosome}", f"{chromosome}")
        output.write(line)

# Updated VCF needs to be bgzipped and tabixed
os.system(f"/root/miniconda2/bin/bgzip {tempFile}")
os.system(f"/root/miniconda2/bin/tabix {tempFile}.gz")
# Phase with Beagle
os.system(f"java -Xmx40g -jar /beagle.25Nov19.28d.jar gt={tempFile}.gz \
out={outputFile} \
chrom={chromosome} \
map=/references/1000GP_Phase3/genetic_map_chr{chromosome}_combined_b37_beagle.txt \
ref=/references/1000GP_Phase3/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
impute=false")

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')