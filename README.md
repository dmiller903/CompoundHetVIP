## Compound Heterozygous Variant Identification Pipeline (CompoundHetVIP)

### Quick Start
CompoundHetVIP is designed to be used with gVCF or VCF files. Please see [ch_pipeline_example.pdf](https://github.com/dmiller903/CompoundHetVIP/blob/master/ch_pipeline_example.pdf) for example code and a description of each step of the pipeline. Here is a brief overview of each step:

1. Keep variant-only sites of VCF or gVCF files

2. Combine each trio into a single file

3. Liftover trio files and individual files from GRCh38 to GRCh37

4. Remove unplaced sites, multiallelic sites, and duplicate sites from lifted files

5. Separate VCF file into chromosome files, then generate plink files for each chromosome file

6. Phase each of the trios with a haplotype reference panel using SHAPEIT2, Beagle, or Eagle2

7. Revert REF/ALT to be congruent with reference panel

8. Concat and merge phased trio chromosome files into one VCF file

9. Trim and normalize VCF file

10. Annotate with snpEff

11. Load VCF as GEMINI database

12. Query for *CH* variants

13. Add Gene Damage Index Scores and Gene lengths to files

The Docker image, [compound-het-vip](https://hub.docker.com/r/dmill903/compound-het-vip), contains all the tools needed to identify compound heterozygous variants using VCF or gVCF files. Tools available and used in the container include: *Plink2* (1, 2), *Picard* (3), *GATK4* (4), *SAMtools* (5), *BCFtools* (5), *SHAPEIT2* (6), *Beagle* (7), *Eagle2* (8), *vt* (9), *SnpEff* (10), *GEMINI* (11), Gene Damage Index (12), and any necessary dependencies.

### Background
A compound heterozygous variant occurs when a person inherits a variant from one parent within a specific gene and also inherits another variant from the other parent at a different position within the same gene (13). The effect of compound heterozygotic inheritance results in two recessive alleles that may cause disease. To detect these types of variants it is necessary to differentiate between paternally and maternally derived nucleotides. If sequencing has already taken place, computational algorithms can be used to help determine which nucleotides were inherited from each parent through a process termed “phasing” (14). 

Phasing requires specific file types which may vary depending on phasing software. Many phasing programs require that input files have been aligned to a specific reference genome, do not contain multiallelic positions, are free of repeat positions, and that each chromosome is phased separately. Figuring out how to prepare files for phasing can be challenging as passing files from program to program may invoke unforeseen incompatibilities. Also, installing specific programs can be challenging because many programs require various dependencies.

### Our Methodology
We have designed our Compound Heterozygous Variant Identification Pipeline (CompoundHetVIP) to overcome many of the time-consuming challenges that researchers may face when trying to phase patient data. By providing reproducible scripts and a Docker Image where these scripts are executed, other researchers will be able to examine our methodology in detail and apply it to their data. Encapsulating our code within containers helps control what software versions are used, what system libraries are used, and creates a cohesive computational environment.

### Contact
If you encounter an issue, please add to the [issue page](https://github.com/dmiller903/CompoundHetVIP/issues)

### References
1. 	S. Purcell, C. Chang, PLINK 2.0 (www.cog-genomics.org/plink/2.0/).
2. 	C. C. Chang, C. C. Chow, L. C. Tellier, S. Vattikuti, S. M. Purcell, J. J. Lee, Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 4, 7 (2015).
3. 	Picard Tools, (available at http://broadinstitute.github.io/picard/).
4. 	R. Poplin, V. Ruano-Rubio, M. A. DePristo, T. J. Fennell, M. O. Carneiro, G. A. Van der Auwera, D. E. Kling, L. D. Gauthier, A. Levy-Moonshine, D. Roazen, K. Shakir, J. Thibault, S. Chandran, C. Whelan, M. Lek, S. Gabriel, M. J. Daly, B. Neale, D. G. MacArthur, E. Banks, Scaling accurate genetic variant discovery to tens of thousands of samples. bioRxiv (2017), p. 201178.
5. 	H. Li, A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 27, 2987–2993 (2011).
6. 	O. Delaneau, B. Howie, A. J. Cox, J.-F. Zagury, J. Marchini, Haplotype estimation using sequencing reads. Am. J. Hum. Genet. 93, 687–696 (2013).
7. 	S. R. Browning, B. L. Browning, Rapid and accurate haplotype phasing and missing-data inference for whole-genome association studies by use of localized haplotype clustering. Am. J. Hum. Genet. 81, 1084–1097 (2007).
8. 	P.-R. Loh, P. Danecek, P. F. Palamara, C. Fuchsberger, Y. A Reshef, H. K Finucane, S. Schoenherr, L. Forer, S. McCarthy, G. R. Abecasis, R. Durbin, A. L Price, Reference-based phasing using the Haplotype Reference Consortium panel. Nat. Genet. 48, 1443–1448 (2016).
9. 	A. Tan, G. R. Abecasis, H. M. Kang, Unified representation of genetic variants. Bioinformatics. 31, 2202–2204 (2015).
10. P. Cingolani, A. Platts, L. L. Wang, M. Coon, T. Nguyen, L. Wang, S. J. Land, X. Lu, D. M. Ruden, A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly . 6, 80–92 (2012).
11. U. Paila, B. A. Chapman, R. Kirchner, A. R. Quinlan, GEMINI: integrative exploration of genetic variation and genome annotations. PLoS Comput. Biol. 9, e1003153 (2013).
12. Y. Itan, L. Shang, B. Boisson, E. Patin, A. Bolze, M. Moncada-Vélez, E. Scott, M. J. Ciancanelli, F. G. Lafaille, J. G. Markle, R. Martinez-Barricarte, S. J. de Jong, X.-F. Kong, P. Nitschke, A. Belkadi, J. Bustamante, A. Puel, S. Boisson-Dupuis, P. D. Stenson, J. G. Gleeson, D. N. Cooper, L. Quintana-Murci, J.-M. Claverie, S.-Y. Zhang, L. Abel, J.-L. Casanova, The human gene damage index as a gene-level approach to prioritizing exome variants. Proc. Natl. Acad. Sci. U. S. A. 112, 13615–13620 (2015).
13. T. Kamphans, P. Sabri, N. Zhu, V. Heinrich, S. Mundlos, P. N. Robinson, D. Parkhomchuk, P. M. Krawitz, Filtering for compound heterozygous sequence variants in non-consanguineous pedigrees. PLoS One. 8, e70151 (2013).
14. Y. Choi, A. P. Chan, E. Kirkness, A. Telenti, N. J. Schork, Comparison of phasing strategies for whole human genomes. PLoS Genet. 14, e1007308 (2018).
