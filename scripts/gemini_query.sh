#!/bin/bash

gemini comp_hets --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_ch_high_med_aaf_cadd.tsv


gemini de_novo --columns "chrom, gene, start, end, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = 'HIGH' or is_lof = 1) or (impact_severity = 'MED' and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated_cadd.db \
> idiopathic_scoliosis/gVCF/idiopathic_scoliosis_de_novo_high_med_aaf_cadd.tsv