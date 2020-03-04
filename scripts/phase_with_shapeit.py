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
parser = argparse.ArgumentParser(description="The parameters for phasing are set so that SHAPEIT2 uses family \
relationship genotype information and also uses a haplotype reference panel")

parser.add_argument('input_file', help='Input file. Use Plink files (bed, bim, fam) if trio, VCF if not trio \
(do not include suffix if plink files)')
parser.add_argument('output_file', help='Name of output file (without suffix)')
parser.add_argument('chromosome_number', help='SHAPEIT2 phases chromosome by chromosome, so it needs to be known which \
chromosome to phase')
parser.add_argument('--is_trio', help='If the VCF files are not trios, indicate with "n".', default='y')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_file
outputFile = args.output_file
chromosome = 'chr' + args.chromosome_number
trio = args.is_trio

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

# Phase using family information and haplotype reference panel
if trio == "y":
    # Check for alignment issues between sample and reference panel
    os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check -B {inputFile} \
    --output-log {outputFile}_check \
    -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
    --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz \
    /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
    /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3")
    # If alignment issues are found, remove problematic positions while phasing using the exclude file from check step
    if os.path.exists(f"{outputFile}_check.snp.strand.exclude"):
        os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {inputFile} --output-log {outputFile}.log \
        -O {outputFile} \
        -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz \
        /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample \
        --thread 3 --no-mcmc --exclude-snp {outputFile}_check.snp.strand.exclude \
        --force --seed 123456789")
    else:
        os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {inputFile} --output-log {outputFile}.log \
        -O {outputFile} \
        -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz \
        /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc \
        --force --seed 123456789")

# Phase only using haplotype reference panel
elif trio == "n":
    # Check for alignment issues between sample and reference panel
    os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check -V {inputFile} \
    --output-log {outputFile}_check \
    -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
    --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz \
    /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
    /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3")
    # If alignment issues are found, remove problematic positions while phasing using the exclude file from check step
    if os.path.exists(f"{outputFile}_check.snp.strand.exclude"):
        os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -V {inputFile} --output-log {outputFile}.log \
        -O {outputFile} \
        -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz \
        /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc \
        --exclude-snp {outputFile}_check.snp.strand.exclude \
        --force --seed 123456789")
    # If no alignment issues are found, do not remove positions while phasing
    else:
        os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -V {inputFile} --output-log {outputFile}.log \
        -O {outputFile} \
        -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz \
        /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc \
        --force --seed 123456789")

# Convert phased files to vcf files and bgzip output vcf
os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert --input-haps {outputFile} \
--output-log {outputFile}_vcf.log --output-vcf {outputFile}.vcf")
os.system(f"bgzip -f {outputFile}.vcf")

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')