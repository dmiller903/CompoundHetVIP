import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Query data for compound heterozygous and de novo variants using various \
filters.")

parser.add_argument('input_database', help='GEMINI database')
parser.add_argument('output_file_prefix', help='Prefix of output file, including path')
parser.add_argument('--num_cores', help='Up to six cores can be used to run the queries in parallel', default=1)

args = parser.parse_args()

#Create variables of each argument from argparse
inputDatabase = args.input_database
outputName = args.output_file_prefix
numCores = int(args.num_cores)
med = "'MED'"
high = "'HIGH'"

# Create function to run queries in parallel
def queries(query):
    os.system(query)

# Create query variables and add them to a list
queryList = []

# Filter CH Variants with minor allele frequency <= 0.005 and cadd >= 20 (Stringent)
query1 = f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_ch_impactHM_aaf005_cadd20.tsv'
queryList.append(query1)

# Filter CH Variants with minor allele frequency <= 0.01 and cadd >= 15 (Less Stringent)
query2 = f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_ch_impactHM_aaf01_cadd15.tsv'
queryList.append(query2)

# Filter CH Variants with no regard for minor allele frequency and cadd >= 15 (No MAF)
query3 = f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_ch_impactHM_cadd15.tsv'
queryList.append(query3)

# Filter for de novo variants with minor allele frequncey <= 0.005 and cadd >= 20 (Stringent)
query4 = f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_aaf005_cadd20.tsv'
queryList.append(query4)

# Filter for de novo variants with minor allele frequncey <= 0.01 and cadd >= 15 (Less Stringent)
query5 = f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_aaf01_cadd15.tsv'
queryList.append(query5)

# Filter for de novo variants with no regard for minor allele frequency and cadd >= 15 (No MAF)
query6 = f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_cadd15.tsv'
queryList.append(query6)

# Use queryList to run all queries in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    executor.map(queries, queryList)

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')