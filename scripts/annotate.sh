#!/bin/bash
#java -Xmx40g -jar /snpEff/snpEff.jar GRCh37.75 -v \
	#idiopathic_scoliosis/gVCF/idopathic_scoliosis_phased_samples_vt.vcf > \
	#idiopathic_scoliosis/gVCF/idopathic_scoliosis_phased_samples_annotated.vcf

java -Xmx40g -jar /snpEff/SnpSift.jar dbNSFP -v \
	-db /snpEff/./data/dbNSFP2.9.txt.gz \
	idiopathic_scoliosis/gVCF/idopathic_scoliosis_phased_samples_annotated.vcf > \
	idiopathic_scoliosis/gVCF/idopathic_scoliosis_phased_samples_annotated_dbNSFP.vcf
