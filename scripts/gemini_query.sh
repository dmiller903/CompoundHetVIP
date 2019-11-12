#!/bin/bash

gemini comp_hets --columns "chrom, gene, start, end, impact_severity" \
--filter "impact_severity = 'HIGH'" \
idiopathic_scoliosis/gVCF/idopathic_scoliosis_phased_samples_annotated.db \
> idiopathic_scoliosis/gVCF/idopathic_scoliosis_ch_high.tsv


#gemini comp_hets --columns "chrom, gene, start, end, impact_severity" \
#--filter "impact_severity = 'HIGH' or impact_severity = 'MED'" \
#idiopathic_scoliosis/gVCF/idopathic_scoliosis_phased_samples_annotated.db \
#> idiopathic_scoliosis/gVCF/idopathic_scoliosis_ch_high_med.tsv