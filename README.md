## Compound Heterozygous Variant Identification Pipeline (CompHetVIP)

### Background
A compound heterozygous variant occurs when a person inherits a variant from one parent within a specific gene and also inherits another variant from the other parent at a different position within the same gene. The effect of compound heterozygotic inheritance results in two recessive alleles that may cause disease. To detect these types of variants it is necessary to differentiate between paternally and maternally derived nucleotides. If sequencing has already taken place, computational algorithms can be used to help determine which nucleotides were inherited from each parent through a process termed “phasing”. 

Phasing requires specific file types which may vary depending on phasing software. Many phasing programs require that input files have been aligned to a specific reference genome, do not contain multiallelic positions, are free of repeat positions, and that each chromosome is phased separately. Figuring out how to prepare files for phasing can be challenging as passing files from program to program may invoke unforeseen incompatibilities. Also, installing specific programs can be challenging because many programs require various dependencies.

### Our Methodology
We have designed our Compound Heterozygous Variant Pipeline to overcome many of the time-consuming challenges that researchers may face when trying to phase patient data. By incorporating our computational pipeline into Docker images and providing reproducible scripts of our analyses, other researchers will be able to examine our methodology in detail and apply it to their data. Encapsulating our code within containers will help control what software versions are used, what system libraries are used, and create a cohesive computational environment.

This pipeline is designed to be used with gVCF or VCF files. Please see [ch_pipeline_example.pdf](https://github.com/dmiller903/ch-pipeline/blob/master/ch_pipeline_example.pdf) for example code and a description of each step of the pipeline.
