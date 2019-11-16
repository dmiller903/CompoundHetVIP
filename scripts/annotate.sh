#!/bin/bash
java -Xmx40g -jar /snpEff/snpEff.jar GRCh37.75 -v \
	idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_vt.vcf > \
	idiopathic_scoliosis/gVCF/idiopathic_scoliosis_phased_samples_annotated.vcf