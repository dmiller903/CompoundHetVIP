import glob
import re
import os
import concurrent.futures
import subprocess

files = []
for file in glob.glob("FM_YTPYPGH1_test/FM_YTPYPGH1_chr*.vcf"):
    matches = re.findall(r"([\w\-\/]+)\.vcf", file)
    files.append(matches[0])

with open("phase_with_shapeit_individual_trio_check.log", "w") as outputFile:
    #def runShapeit(file):
    for file in files:
        chr = re.findall(r"[\w\-]+\/[\w\-]+(chr[0-9][0-9]?)", file)[0]
        outputProcess = subprocess.Popen("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check -B {} --output-log {}.checks -M /references/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_{}_combined_b37.txt \
        --input-ref /references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3.sample /references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3_{}_impute_macGT1.hap.gz \
        /references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3_{}_impute_macGT1.legend.gz".format(file, file, chr, chr, chr), shell = True, stdout = subprocess.PIPE)
        stdOut, stdErr = outputProcess.communicate()
        stdOut = stdOut.decode("utf-8")
        outputFile.write(stdOut)
        stdErr = stdErr.decode("utf-8")
        outputFile.write(stdErr)
        #os.system("/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {} -M /references/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_{}_combined_b37.txt \
                #--input-ref /references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3.sample /references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3_{}_impute_macGT1.hap.gz \
                #/references/ALL_1000G_phase1integrated_v3_impute_macGT1/ALL_1000G_phase1integrated_v3_{}_impute_macGT1.legend.gz -O {}_phased --output-graph {}_phased.graph \
                #--output-log {}_phased.log --exclude-snp /fslhome/dmill903/research/idiopathic_scoliosis/gVCF/FM_YTPYPGH1_test/FM_YTPYPGH1_chr21.checks.snp.strand.exclude".format(file, chr, chr, chr, file, file, file))
        #os.system("shapeit -convert --input-haps {}.phased --output-vcf {}.phased.vcf".format(file, file))

    #with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        #executor.map(runShapeit, files)
        #executor.close()
