import re
import os
import glob
import concurrent.futures

#Update original genetic map files format to Eagle and Beagle formats
def updateFiles(file):
    fileName = re.findall(r"([\w/_]+genetic_map_chr(\w+)_combined_b37)\.txt", file)[0][0]
    chrom = re.findall(r"([\w/_]+genetic_map_chr(\w+)_combined_b37)\.txt", file)[0][1]
    eagleOutput = f"{fileName}_eagle.txt"
    beagleOutput = f"{fileName}_beagle.txt"
    with open(file) as inputFile, open(eagleOutput, 'w') as eagleOutput, open(beagleOutput, 'w') as beagleOutput:
        header = inputFile.readline()
        header = "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
        eagleOutput.write(header)
        for line in inputFile:
            #Eagle format
            eagleLine = f"{chrom} " + line
            eagleOutput.write(eagleLine)
            #Beagle format
            beagleLineList = line.rstrip().split(" ")
            beagleLine = f"{chrom} {beagleLineList[1]} {beagleLineList[2]} {beagleLineList[0]}\n"
            beagleOutput.write(beagleLine)
    #Download reference VCF's
    os.system(f"wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    && mv ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz /references/1000GP_Phase3/")
    os.system(f"wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi \
    && mv ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi /references/1000GP_Phase3/")

with concurrent.futures.ProcessPoolExecutor(max_workers=22) as executor:
    executor.map(updateFiles, glob.glob("/references/1000GP_Phase3/genetic_map_chr*_combined_b37.txt"))