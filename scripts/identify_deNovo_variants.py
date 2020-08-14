import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Uses a GEMINI database as input to identify de Novo variants. If a \
--fam_file' is not used, false positives may result since parental haplotypes not being takin in to consideration.")

parser.add_argument('input_file', help='GEMINI database')
parser.add_argument('output_file', help='Name of output file')
parser.add_argument('fam_file', help='Family relationships are mandatory to help with the de Novo identification process.')
parser.add_argument('--cadd', help='If you use strict argument, and want to customize cadd cut-off value.', default='15')
parser.add_argument('--maf', help='If you use strict argument, and want to customize maf cut-off value.', default='0.01')

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_file
outputFile = args.output_file
familyFile = args.fam_file
inputCadd = float(args.cadd)
inputMaf = args.maf
if inputMaf != "None":
    inputMaf = float(args.maf)

#Function to get convert sample genotype from alpha to numeric
def getNumericGenotype(genotype, ref, alt):
    if "|" in genotype and "." not in genotype:
        genotypeList = genotype.split("|")
        firstAllele = ""
        secondAllele = ""
        if genotypeList[0] == ref:
            firstAllele = "0"
        elif genotypeList[0] == alt:
            firstAllele = "1"
        else:
            firstAllele = "."
        if genotypeList[1] == ref:
            secondAllele = "0"
        elif genotypeList[1] == alt:
            secondAllele = "1"
        else:
            secondAllele = "."
        newGenotype = f"{firstAllele}|{secondAllele}"
        return(newGenotype)
    else:
        return(".|.")

#Function to grab information from header of input file
def getHeaderInfo(headerList):
    startIndex = headerList.index("start")
    geneIndex = headerList.index("gene")
    refIndex = headerList.index("ref")
    altIndex = headerList.index("alt")
    impactIndex = headerList.index("impact_severity")
    caddIndex = headerList.index("cadd_scaled")
    mafIndex = headerList.index("aaf_1kg_all")
    lofIndex = headerList.index("is_lof")
    exonicIndex = headerList.index("is_exonic")
    samples = headerList[13:]
    return(startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, mafIndex, lofIndex, exonicIndex, samples)

#Function to grab information from line of input file
def getLineInfo(lineList):
    start = lineList[startIndex]
    gene = lineList[geneIndex]
    ref = lineList[refIndex]
    alt = lineList[altIndex]
    impact = lineList[impactIndex]
    cadd = lineList[caddIndex]
    maf = lineList[mafIndex]
    lof = lineList[lofIndex]
    exonic = lineList[exonicIndex]
    return(start, gene, ref, alt, impact, cadd, maf, lof, exonic)

def iterateThroughSamples():
    for sampleIndex in sampleIndexes:
        sample = headerList[sampleIndex]
        genotype = lineList[sampleIndex]
        newGenotype = getNumericGenotype(genotype, ref, alt)
        if gene not in sampleGenotype[sample] and "." not in newGenotype:
            sampleGenotype[sample][gene] = [newGenotype]
            samplePositions[sample][gene] = [start]
        elif gene in sampleGenotype[sample] and "." not in newGenotype:
            sampleGenotype[sample][gene].append(newGenotype)
            samplePositions[sample][gene].append(start)

# Create a .tsv that has all pertinent information for compound heterozygous identification
impactSeverity = "'LOW'"
geminiTsv = f"{inputFile.replace('.db', '_gemini.tsv')}"
if not os.path.exists(geminiTsv):
    os.system(f'gemini query --header -q "select chrom, start, vcf_id, ref, alt, gene, is_exonic, impact_severity, \
        is_lof, aaf_1kg_all, cadd_scaled, impact, biotype, (gts).(*) from variants where impact_severity != {impactSeverity}" \
        {inputFile} \
        > {geminiTsv}')

# Use fam file to create a list of samples, list of parents, and a parent dictionary where each key is a parent ID and value is sample ID
parentDict = {}
parentList = []
patientList = []
familyDict = {}
with open(familyFile) as familyFile:
    for line in familyFile:
        lineList = line.rstrip("\n").split("\t")
        if lineList[-1] is "2":
            parentDict[f"gts.{lineList[2]}"] = f"gts.{lineList[1]}"
            parentDict[f"gts.{lineList[3]}"] = f"gts.{lineList[1]}"
            patientList.append(f"gts.{lineList[1]}")
            familyDict[f"gts.{lineList[1]}"] = [f"gts.{lineList[2]}", f"gts.{lineList[3]}"]
        else:
            parentList.append(f"gts.{lineList[1]}")

"""
Iterate through inputFile in order to create two dictionaries: sampleGenotype and samplePositions. The key of 
sampleGenotype is the sample ID and value is a dictionary where the key is a gene and the value is a list of all 
genotypes ("0|1", "1|0", or "0|0") for that gene that meet specific CADD score, minor allele frequency, and impact 
severity criteria. The samplePositions has the same information, except the list for each gene is genotype positions.
"""

sampleGenotype = {}
samplePositions = {}
sampleIndexes = []
with open(geminiTsv) as geminiFile:
    header = geminiFile.readline()
    headerList = header.rstrip("\n").split("\t")
    startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, mafIndex, lofIndex, exonicIndex, samples = getHeaderInfo(headerList)
    for sample in samples:
        sampleIndexes.append(headerList.index(sample))
        sampleGenotype[sample] = {}
        samplePositions[sample] = {}    
    for line in geminiFile:
        lineList = line.rstrip("\n").split("\t")
        start, gene, ref, alt, impact, cadd, maf, lof, exonic = getLineInfo(lineList)
        if cadd != "None" and maf != "None":
            if ((impact == "HIGH" or lof == "1") or (impact == "MED" and float(cadd) >= inputCadd)) and float(maf) <= inputMaf:
                iterateThroughSamples()
        elif cadd == "None" and maf == "None":
            if impact == "HIGH" or lof == "1":
                iterateThroughSamples()
        elif cadd != "None" and maf == "None":
            if (impact == "HIGH" or lof == "1") or (impact == "MED" and float(cadd) >= inputCadd):
                iterateThroughSamples()
        elif cadd == "None" and maf != "None":
            if (impact == "HIGH" or lof == "1") and float(maf) <= inputMaf:
                iterateThroughSamples()
print("Sample Dictionaries Created.")

"""
Use sampleGenotype and samplePositions to generate a new dictionaries where the key is the sample ID and the value is 
a dictionary where the key is a gene and the value is a list of genotypes (or positions) where de Novo variant(s) are found.
"""
deNovoPositionDict = {}
deNovoGenotypeDict = {}
for patient in patientList:
    deNovoPositionDict[patient] = {}
    deNovoGenotypeDict[patient] = {}
    parent1 = familyDict[patient][0]
    parent2 = familyDict[patient][1]
    for gene, genotypes in sampleGenotype[patient].items():
        for i, genotype in enumerate(genotypes):
            positionList = samplePositions[patient][gene]
            position = positionList[i]
            #This part helps eliminate genotypes being added to the deNovo list where either parent is homozygous recessive
            parentGenotype1 = ""
            parentGenotype2 = ""
            if gene in samplePositions[parent1] and position in samplePositions[parent1][gene]:
                parentPosIndex = samplePositions[parent1][gene].index(position)
                parentGenotype1 = sampleGenotype[parent1][gene][parentPosIndex]
            if gene in samplePositions[parent2] and position in samplePositions[parent2][gene]:
                parentPosIndex = samplePositions[parent2][gene].index(position)
                parentGenotype2 = sampleGenotype[parent2][gene][parentPosIndex]
            if "1" not in parentGenotype1 and "1" not in parentGenotype2:
                if genotype in ["1|0", "0|1", "1|1"] and gene not in deNovoPositionDict[patient]:
                    deNovoPositionDict[patient][gene] = [position]
                    deNovoGenotypeDict[patient][gene] = [genotype]
                elif genotype in ["1|0", "0|1", "1|1"] and gene in deNovoPositionDict[patient]:
                    deNovoPositionDict[patient][gene].append(position)
                    deNovoGenotypeDict[patient][gene].append(genotype)
print("de Novo variant dictionaries created.")

#Iterate through the input file and use the deNovoPositionDict in order to output de Novo variant data for each sample
with open(geminiTsv) as geminiFile, open(outputFile, "w") as outputFile:
    header = geminiFile.readline()
    headerList = header.rstrip("\n").split("\t")
    startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, mafIndex, lofIndex, exonicIndex, samples = getHeaderInfo(headerList)
    columnInfo = headerList[0:13]
    newHeader = "\t".join(columnInfo) + "\tgenotype\tsample\n"
    outputFile.write(newHeader)
    sampleIndexes = []
    for patient in patientList:
        patientIndex = headerList.index(patient)
        sampleIndexes.append(patientIndex)
    for line in geminiFile:
        lineList = line.rstrip("\n").split("\t")
        start, gene, ref, alt, impact, cadd, maf, lof, exonic = getLineInfo(lineList)
        for sampleIndex in sampleIndexes:
            sample = headerList[sampleIndex]
            parent1 = familyDict[sample][0]
            parent2 = familyDict[sample][1]
            if gene in deNovoPositionDict[sample] and start in deNovoPositionDict[sample][gene]:
                genotype = lineList[sampleIndex]
                numericGenotype = getNumericGenotype(genotype, ref, alt)
                if "." not in numericGenotype and numericGenotype in ["1|0", "0|1", "1|1"]:
                    columnInfo = lineList[0:13]
                    columnStr = "\t".join(columnInfo)
                    newLine = f"{columnStr}\t{numericGenotype}\t{sample.replace('gts.', '')}\n"
                    outputFile.write(newLine)

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')