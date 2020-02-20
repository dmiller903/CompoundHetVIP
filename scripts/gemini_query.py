import re
import os
from sys import argv

#Get disease name based on path
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

diseaseName = chromosome = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]
med = "'MED'"
high = "'HIGH'"

# Filter CH Variants with minor allele frequency <= 0.005 and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf005_cadd20.tsv')

# Filter CH Variants with minor allele frequency <= 0.01 and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf01_cadd20.tsv')

# Filter CH Variants with no regard for minor allele frequency and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_cadd20.tsv')

# Filter CH Variants with minor allele frequency <= 0.005 and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf005_cadd15.tsv')

# Filter CH Variants with minor allele frequency <= 0.01 and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf01_cadd15.tsv')

# Filter CH Variants with no regard for minor allele frequency and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_cadd15.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.005 and cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_aaf005_cadd20.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.01 and cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_aaf01_cadd20.tsv')

# Filter for de novo variants with no regard for minor allele frequency and cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_cadd20.tsv')