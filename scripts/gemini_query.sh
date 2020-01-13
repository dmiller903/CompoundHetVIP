#!/bin/bash

# Filter CH Variants with minor allele frequency <= 0.005
gemini comp_hets --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_ch_impactHM_aaf005_cadd20.tsv

# Filter CH Variants with minor allele frequency <= 0.01
gemini comp_hets --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_ch_impactHM_aaf01_cadd20.tsv

# Filter CH Variants with no regard for minor allele frequency
gemini comp_hets --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_ch_impactHM_cadd20.tsv

# Filter CH Variants with CADD score >= 20
gemini comp_hets --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "cadd_scaled >=20" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_ch_cadd20.tsv

# Filter CH Variants with no filters
gemini comp_hets --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_ch_noFilters.tsv

# Filter for de novo variants with minor allele frequncey <= 0.005
gemini de_novo --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_de_novo_impactHM_aaf005_cadd20.tsv

# Filter for de novo variants with minor allele frequncey <= 0.01
gemini de_novo --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_de_novo_impactHM_aaf01_cadd20.tsv

# Filter for de novo variants with no regard for minor allele frequency
gemini de_novo --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_de_novo_impactHM_cadd20.tsv