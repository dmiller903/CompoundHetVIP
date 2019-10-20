import glob
import re
import os
import concurrent.futures

folders = glob.glob("combined_trio/separate_by_chr/chr[0-9][0-9]") + glob.glob("combined_trio/separate_by_chr/chr[0-9]")
files = []
for folder in folders:
    for file in glob.glob("{}/*.vcf".format(folder)):
        matches = re.findall("(.+)\.vcf", file)
        files.append(matches[0])

def runShapeit(file):
    chr = re.findall(".+/(chr[0-9][0-9]?)", file)[0]
    folder = re.findall("(.+/)chr[0-9][0-9]?", file)[0]
    #os.system("shapeit -check -B {} --output-log {}.checks".format(file, file))
    os.system("shapeit -B {} -M /references/1000GP_Phase3/genetic_map_{}_combined_b37.txt -O {}_family_phased \
    --output-graph {}_family.graph --output-log {}_phased.log".format(file, chr, file, file, file))
    os.system("shapeit -convert --input-haps {}.phased --output-vcf {}.phased.vcf".format(file, file))

with concurrent.futures.ProcessPoolExecutor(max_workers=22) as executor:
    executor.map(runShapeit, files)

