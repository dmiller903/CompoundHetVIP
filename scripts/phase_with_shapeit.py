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
parser = argparse.ArgumentParser(description="The parameters for phasing are set so that SHAPEIT2 uses family \
relationship genotype information and also uses a haplotype reference panel")

parser.add_argument('input_plink', help='Input file (without suffix)')
parser.add_argument('output_file', help='Name of output file (without suffix)')
parser.add_argument('chromosome_number', help='Shapeit phases chromosome by chromsome, so it needs to be known which \
chromosome to phase')
parser.add_argument('--is_trio', help='If the VCF files are not trios, indicate with "n".', default='y')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_plink
outputFile = args.output_file
chromosome = 'chr' + args.chromosome_number
trio = args.is_trio

if trio == "y":
    os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check -B {} --output-log {}_check \
    -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
    --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
    /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3".format(inputFile, outputFile, chromosome, chromosome, chromosome))
    if os.path.exists("{}_check.snp.strand.exclude".format(outputFile)):
        os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {} --output-log {}.log -O {} \
        -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc --exclude-snp {}_check.snp.strand.exclude \
        --force --seed 123456789".format(inputFile, outputFile, outputFile, chromosome, chromosome, chromosome, outputFile))
    else:
        os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {} --output-log {}.log -O {} \
        -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc \
        --force --seed 123456789".format(inputFile, outputFile, outputFile, chromosome, chromosome, chromosome))

elif trio == "n":
    os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check -V {} --output-log {}_check \
    -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
    --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
    /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3".format(inputFile, outputFile, chromosome, chromosome, chromosome))
    if os.path.exists("{}_check.snp.strand.exclude".format(outputFile)):
        os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -V {} --output-log {}.log -O {} \
        -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc --exclude-snp {}_check.snp.strand.exclude \
        --force --seed 123456789".format(inputFile, outputFile, outputFile, chromosome, chromosome, chromosome, outputFile))
    else:
        os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -V {} --output-log {}.log -O {} \
        -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc \
        --force --seed 123456789".format(inputFile, outputFile, outputFile, chromosome, chromosome, chromosome))
    
os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert --input-haps {} \
--output-log {}_vcf.log --output-vcf {}.vcf".format(outputFile, outputFile, outputFile))
os.system(f"bgzip -f {outputFile}.vcf")