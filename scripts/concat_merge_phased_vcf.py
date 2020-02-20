import os
import time
import argparse
import re
import gzip
import glob
#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Phased results can have the REF and ALT alleles switched as compared \
to the reference genome. We are unsure exactly why this occurs. For files with trios, it may be that since each file \
only has 3 samples, the ALT allele is more common in the trio and becomes the REF. This step ensures that the REF/ALT \
alleles of the phased VCF files are congruent with the REF/ALT of the reference genome. In addition, sites with Mendel \
errors are removed.")

parser.add_argument('input_vcf', help='Input file')
parser.add_argument('output_file', help='Name of output file')
parser.add_argument('chromosome_number', help='Chromosome number is needed so the script can determine which reference \
file to use.')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file
chromosome = args.chromosome_number
tempFile = "/tmp/" + re.findall(r'/?([\w\-_\.]+)', outputFile)[-1]

def concatMerge(trio):
    files = fileDict[trio]

    for index, file in enumerate(files):
        #os.system("gzip -d {}.gz".format(file))
        os.system("bgzip -f {} && tabix -fp vcf {}.gz".format(file, file))
        files[index] = "{}.gz".format(file)


    fileName = re.findall(r"([\w\-\/_]+\/[\w\-_]+)_chr[A-Z0-9][A-Z0-9]?_phased_reverted\.vcf", files[0])[0]
    outputName = "{}_phased_combined.vcf".format(fileName)
    files = " ".join(files)
    os.system("bcftools concat {} -o {}".format(files, outputName))
    os.system("bgzip -f {} && tabix -fp vcf {}.gz".format(outputName, outputName))

with concurrent.futures.ProcessPoolExecutor(max_workers=35) as executor:
    executor.map(concatMerge, fileDict)

# Merge all phased, concatenated, trio files into one    
concatFilesString = " ".join(concatFiles)
outputName = "{}/{}_phased_samples.vcf".format(pathToFiles, diseaseName)
os.system("bcftools merge -m both {} -o {}".format(concatFilesString, outputName))
os.system("bgzip -f {} && tabix -fp vcf {}.gz".format(outputName, outputName))

# Create a merged family file
# create a proband dictionary where the key is the sampleId and the value is the familyId
# also create a parent dictionary where the key is familyId and the value is a dictionary that has a key of the sampleId and value of gender
probandDict = {}
parentDict = {}
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    probandIndex = headerList.index("proband")
    genderIndex = headerList.index("sex")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        probandStatus = sampleData[probandIndex]
        gender = sampleData[genderIndex]
        if probandStatus == "Yes":
            probandDict[sampleId] = sampleFamilyId
        else:
            if sampleFamilyId not in parentDict:
                parentDict[sampleFamilyId] = {sampleId: gender}
            else:
                parentDict[sampleFamilyId][sampleId] = gender

# Create a dictionary where each sample has the rest of the family information needed for the family file
sampleDict = dict()
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    probandIndex = headerList.index("proband")
    genderIndex = headerList.index("sex")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        probandStatus = sampleData[probandIndex]
        gender = sampleData[genderIndex]
        paternal = ""
        maternal = ""
        if probandStatus == "Yes":
            familyDict = parentDict[sampleFamilyId]
            for key, value in familyDict.items():
                if value == "1":
                    paternal = key
                else:
                    maternal = key
            sampleDict[sampleId] = "{}\t{}\t{}\t{}\t{}\t2\n".format(sampleFamilyId, sampleId, paternal, maternal, gender)
        else:
            sampleDict[sampleId] = "{}\t{}\t0\t0\t{}\t1\n".format(sampleFamilyId, sampleId, gender)
            
# create a sample list in the order of the vcf file
with gzip.open("{}/{}_phased_samples.vcf.gz".format(pathToFiles, diseaseName), "rt") as vcfFile:
    for line in vcfFile:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            sampleList = line.rstrip().split("\t")[9:]
        else:
            break

# use the sample order in the list to output each sample in order as found in the vcf file
with open("{}/{}.fam".format(pathToFiles, diseaseName), "w") as outputFile:
    for sample in sampleList:
        outputFile.write(sampleDict[sample])