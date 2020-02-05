from sys import argv
import re
import statistics

inputFile = argv[1]
outputFile = argv[2]
ncbiFile = argv[3]
ensembleFile = argv[4]

# Create new sample names based on family id
familyDict = {}
with open(inputFile) as chFile:
    header = chFile.readline()
    headerList = header.rstrip("\n").split("\t")
    familyIndex = headerList.index("family_id")
    familyCount = 1
    for line in chFile:
        lineList = line.rstrip("\n").split("\t")
        familyId = lineList[familyIndex]
        if familyId not in familyDict:
            familyDict[familyId] = "patient_{}".format(familyCount)
            familyCount += 1

# Create a ncbi dictionary of gene lengths
ncbiLength = {}
with open(ncbiFile) as ncbiFile:
    header = ncbiFile.readline()
    headerList = header.rstrip("\n").split("\t")
    geneNameIndex = headerList.index("Symbol")
    aliasesIndex = headerList.index("Aliases")
    startIndex = headerList.index("start_position_on_the_genomic_accession")
    endIndex = headerList.index("end_position_on_the_genomic_accession")
    for line in ncbiFile:
        lineList = line.rstrip("\n").split("\t")
        geneName = lineList[geneNameIndex]
        start = lineList[startIndex]
        end = lineList[endIndex]
        aliases = lineList[aliasesIndex]
        aliases = aliases.rstrip("\n").split(", ")
        if start == "" or end == "":
            ncbiLength[geneName] = [aliases, ["NA"]]
        else:
            length = int(end) - int(start)
            ncbiLength[geneName] = [aliases, [length]]

ensembleLength = {}
# Create an ensemble dictionary of gene lengths
with open(ensembleFile) as ensembleFile:
    for line in ensembleFile:
        if "#" not in line:
            lineList = line.rstrip("\n").split("\t")
            geneName = lineList[8]
            if re.findall(r'gene_name "([\w\.\-]+)"', geneName):
                geneName = re.findall(r'gene_name "([\w\.\-]+)"', geneName)[0]
            else:
                continue
            start = lineList[3]
            end = lineList[4]
            if start == "." or end == ".":
                ensembleLength[geneName] = "NA"
            else:
                length = int(end) - int(start)
                if geneName not in ensembleLength:
                    ensembleLength[geneName] = length
                elif geneName in ensembleLength and length > ensembleLength[geneName]:
                    ensembleLength[geneName] = length


with open(inputFile) as chFile, open(outputFile, "w") as tidyFile:
    header = chFile.readline()
    headerList = header.rstrip("\n").split("\t")
    chromIndex = headerList.index("chrom")
    geneIndex = headerList.index("gene")
    startIndex = headerList.index("start")
    endIndex = headerList.index("end")
    impactSeverityIndex = headerList.index("impact_severity")
    aafIndex = headerList.index("aaf_1kg_all")
    caddIndex = headerList.index("cadd_scaled")
    impactIndex = headerList.index("impact")
    biotypeIndex = headerList.index("biotype")
    familyIndex = headerList.index("family_id")
    tidyFile.write("patient\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgender\t{}\t{}\tgene_length\n".format(headerList[geneIndex], 
    headerList[chromIndex], headerList[impactSeverityIndex], headerList[impactIndex], headerList[biotypeIndex], 
    headerList[aafIndex], headerList[caddIndex], headerList[startIndex], headerList[endIndex]))
    for line in chFile:
        lineList = line.rstrip("\n").split("\t")
        patient = lineList[familyIndex]
        patient = familyDict[patient]
        gene = lineList[geneIndex]
        chrom = lineList[chromIndex]
        impactSeverity = lineList[impactSeverityIndex]
        impact = lineList[impactIndex]
        biotype = lineList[biotypeIndex]
        aaf = lineList[aafIndex]
        cadd = lineList[caddIndex]
        gender = re.findall(r"affected;(\w+)\)", line)[0]
        start = lineList[startIndex]
        end = lineList[endIndex]
        if gene in ncbiLength:
            length = ncbiLength[gene][1][0]
            tidyFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(patient, gene, chrom, impactSeverity, 
            impact, biotype, aaf, cadd, gender, start, end, length))
        elif gene in ensembleLength:
            length = ensembleLength[gene]
            tidyFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(patient, gene, chrom, impactSeverity, 
            impact, biotype, aaf, cadd, gender, start, end, length))
        else:
            inAliases = False
            output = ""
            for key, value in ncbiLength.items():
                if gene in value[0]:
                    length = value[1][0]
                    output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(patient, gene, chrom, impactSeverity, 
                    impact, biotype, aaf, cadd, gender, start, end, length)
                    break
                else:
                    output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA\n".format(patient, gene, chrom, impactSeverity, 
                    impact, biotype, aaf, cadd, gender, start, end)
            tidyFile.write(output)

# Get average gene lengths
ncbiLengthList = list()
for key, value in ncbiLength.items():
    length = value[1][0]
    if isinstance(length, str):
        continue
    ncbiLengthList.append(int(length))
ncbiAverageLength = statistics.mean(ncbiLengthList)
ncbiMedian = statistics.median(ncbiLengthList)
print(ncbiMedian)

ensemblLengthList = list()
for key, value in ensembleLength.items():
    length = value
    if isinstance(length, str):
        continue
    ensemblLengthList.append(int(length))
ensemblAverageLength = statistics.mean(ensemblLengthList)
ensemblMedian = statistics.median(ensemblLengthList)
print(ensemblMedian)

combinedAvg = (ncbiAverageLength + ensemblAverageLength) / 2
combinedMedian = (ncbiMedian + ensemblMedian) / 2
print("The average gene length across the human genome is {}".format(combinedAvg))
print("The median gene length across the human genome is {}".format(combinedMedian))

    
