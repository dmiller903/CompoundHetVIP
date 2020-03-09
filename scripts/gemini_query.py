import os
import time
import argparse

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Annotate VCF using snpEff")

parser.add_argument('input_database', help='GEMINI database')
parser.add_argument('output_file_prefix', help='Prefix of output file, including path')

args = parser.parse_args()

#Create variables of each argument from argparse
inputDatabase = args.input_database
outputName = args.output_file_prefix
med = "'MED'"
high = "'HIGH'"

# CH Variants with no filters
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
{inputDatabase} \
> {outputName}_ch_no_filter.tsv')

# Filter CH Variants with minor allele frequency <= 0.005 and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_ch_impactHM_aaf005_cadd20.tsv')

# Filter CH Variants with minor allele frequency <= 0.01 and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_ch_impactHM_aaf01_cadd20.tsv')

# Filter CH Variants with no regard for minor allele frequency and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_ch_impactHM_cadd20.tsv')

# Filter CH Variants with minor allele frequency <= 0.005 and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_ch_impactHM_aaf005_cadd15.tsv')

# Filter CH Variants with minor allele frequency <= 0.01 and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_ch_impactHM_aaf01_cadd15.tsv')

# Filter CH Variants with no regard for minor allele frequency and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_ch_impactHM_cadd15.tsv')

# de novo variants with no filters
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
{inputDatabase} \
> {outputName}_de_novo_no_filter.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.005 and cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_aaf005_cadd20.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.01 and cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_aaf01_cadd20.tsv')

# Filter for de novo variants with no regard for minor allele frequency and cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=20)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_cadd20.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.005 and cadd >= 15
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_aaf005_cadd15.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.01 and cadd >= 15
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_aaf01_cadd15.tsv')

# Filter for de novo variants with no regard for minor allele frequency and cadd >= 15
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=15)" \
{inputDatabase} \
> {outputName}_de_novo_impactHM_cadd15.tsv')

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')