# Import necessary modules
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

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_file', help='Name of output VCF file')
parser.add_argument('chromosome_number', help='Chromosome number is needed so the script can determine which reference \
file to use')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_file.rstrip(".gz")
chromosome = args.chromosome_number
tempFile = "/tmp/" + re.findall(r'/?([\w\-_\.]+)', outputFile)[-1]

# Create dictionary from legend file where key is chromosome and value is a set of strings that give site information
posDict = dict()
for file in glob.glob("/references/1000GP_Phase3/*legend.gz"):
    if file == f"/references/1000GP_Phase3/1000GP_Phase3_chr{chromosome}.legend.gz":
        with gzip.open(file, 'rt') as reference:
            header = reference.readline()
            headerList = header.rstrip().split()
            refIndex = headerList.index("a0")
            altIndex = headerList.index("a1")
            posIndex = headerList.index("position")
            idIndex = headerList.index("id")
            for line in reference:
                lineList = line.rstrip().split(" ")
                pos = lineList[posIndex]
                ref = lineList[refIndex]
                alt = lineList[altIndex]
                siteStr = f"{pos} {ref} {alt}"
                if chromosome not in posDict:
                    posDict[chromosome] = {siteStr}
                else:
                    posDict[chromosome].add(siteStr)

mendelErrorCount = 0
fileWithoutSuffix = re.findall(r'([\w\-_/]+)\.', inputFile)[0]
mendelErrorFile = f"{fileWithoutSuffix}.snp.me"
mendelErrorSet = set()
# Create a set of any positions with mendel errors as given by the shapeit2 .snp.me files
if os.path.exists(mendelErrorFile):
    with open(mendelErrorFile) as mendelFile:
        for line in mendelFile:
            lineSplit = line.split("\t")
            mendelError = lineSplit[2]
            pos = lineSplit[1]
            if mendelError == "1":
                mendelErrorSet.add(pos)

# Flip alt and ref, and remove mendel errors if shapeit2 was used to phase
rawCount = 0
flipCount = 0
total = 0
with gzip.open(inputFile, 'rt') as sample, gzip.open(tempFile, 'wb') as output:
    for line in sample:
        if "##" in line:
            output.write(line.encode())
        elif line.startswith("#CHROM"):
            header = line.split("\t")
            chromIndex = header.index("#CHROM")
            posIndex = header.index("POS")
            refIndex = header.index("REF")
            altIndex = header.index("ALT")
            output.write(line.encode())
        else:
            lineList = line.split("\t")
            chrom = lineList[chromIndex]
            pos = lineList[posIndex]
            ref = lineList[refIndex]
            alt = lineList[altIndex]
            rawStr = f"{pos} {ref} {alt}"
            flipStr = f"{pos} {alt} {ref}"
            if rawStr in posDict[chrom] and pos not in mendelErrorSet:
                output.write(line.encode())
                rawCount += 1
                total += 1
            elif flipStr in posDict[chrom] and pos not in mendelErrorSet:
                lineList[refIndex] = alt
                lineList[altIndex] = ref
                line = "\t".join(lineList)
                line = line.replace("0|1", "b|a").replace("1|0", "a|b").replace("1|1", "a|a").replace("0|0", "b|b")
                line = line.replace("b|a", "1|0").replace("a|b", "0|1").replace("a|a", "0|0").replace("b|b", "1|1")
                output.write(line.encode())
                flipCount += 1
                total += 1
            else:
                total += 1
                mendelErrorCount += 1

# Output statistics
rawPercent = (rawCount / total) * 100
flipPercent = (flipCount / total) * 100
totalPercent = ((flipCount + rawCount) / total) * 100
if flipCount == 0 and mendelErrorCount == 0:
    print(f"For {inputFile}, chr{chromosome}, {rawCount} ({rawPercent:.2f}%) of the sites were unchanged. No outputFile was generated.")
else:
    os.system(f"mv {tempFile} {outputFile}.gz")
    print(f"For {inputFile}, chr{chromosome}, {rawCount} ({rawPercent:.2f}%) of the sites were unchanged")
    print(f"For {inputFile}, chr{chromosome}, {flipCount} ({flipPercent:.2f}%) of the sites were switched to match the reference panel")
    print(f"For {inputFile}, chr{chromosome}, {totalPercent:.2f}% of the sites are now congruent with the reference panel\n")
    print(f"For {inputFile}, chr{chromosome}, {mendelErrorCount} sites were removed due to mendel errors\n")

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')