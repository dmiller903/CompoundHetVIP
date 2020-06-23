from sys import argv
import time
import argparse

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Keeps CH variants where the combined aaf for each variant is <= a user \
defined cutoff value.')

parser.add_argument('GEMINI_file', help='.tsv file generated using a GEMINI database')
parser.add_argument('output_file', help='Name of output file')
parser.add_argument('cutoff_value', help='The aaf cutoff value')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.GEMINI_file
outputFile = args.output_file
cutoff = float(args.cutoff_value)

"""
Use the input file to generate a dictionary where the key is a sample ID and the value is a dictionary where the key is a gene and the value is a list of aaf values for that gene
"""
with open(inputFile) as geminiFile:
    header = geminiFile.readline()
    headerList = header.rstrip().split('\t')
    aafIndex = headerList.index('aaf_1kg_all')
    sampleIndex = headerList.index('sample')
    geneIndex = headerList.index('gene')
    
    sampleDict = {}
    for line in geminiFile:
        lineList = line.rstrip().split('\t')
        sample = lineList[sampleIndex]
        gene = lineList[geneIndex]
        aaf = lineList[aafIndex]
        if sample not in sampleDict:
            sampleDict[sample] = {gene: [aaf]}
        elif sample in sampleDict and gene not in sampleDict[sample]:
            sampleDict[sample][gene] = [aaf]
        elif sample in sampleDict and gene in sampleDict[sample]:
            sampleDict[sample][gene].append(aaf)
"""
Create a new dictionary where the key is a sample ID and the value is a dictionary where each key is a gene and value is
the product of the aaf's for that gene
"""
updatedSampleDict = {}
for sample, geneDict in sampleDict.items():
    updatedSampleDict[sample] = {}
    for gene, aafs in geneDict.items():
        updatedSampleDict[sample][gene] = 1
        for aaf in aafs:
            updatedSampleDict[sample][gene] *= float(aaf)

"""
Create a new file where only the genes for each patient have a combined aaf >= the cutoff value are kept
"""
with open(inputFile) as geminiFile, open(outputFile, 'w') as outputFile:
    header = geminiFile.readline()
    headerList = header.rstrip().split('\t')
    aafIndex = headerList.index('aaf_1kg_all')
    sampleIndex = headerList.index('sample')
    geneIndex = headerList.index('gene')
    outputFile.write(f'{header.rstrip()}\tcombined_aaf\n')
    sampleDict = {}
    for line in geminiFile:
        lineList = line.rstrip().split('\t')
        sample = lineList[sampleIndex]
        gene = lineList[geneIndex]
        combinedAaf = updatedSampleDict[sample][gene]
        if combinedAaf <= cutoff:
            outputFile.write(f'{line.rstrip()}\t{combinedAaf}\n')

# Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')