import glob
from sys import argv

manifestPath = argv[1]
biospecimenPath = argv[2]
outputPath = argv[3]

initialDict = {}
finalDict = {}

with open(manifestPath) as manifest, open(biospecimenPath) as biospecimen, open(outputPath, 'w') as output:
    manifestColumnNames = manifest.readline()
    manifestColumnNames = manifestColumnNames.rstrip().split("\t")
    fileNameIndex = manifestColumnNames.index("File Name")
    familyIdIndex = manifestColumnNames.index("Family Id")
    manifestExternalIdIndex = manifestColumnNames.index("Aliquot External ID")
    probandIndex = manifestColumnNames.index("Proband")
    for line in manifest:
        line = line.rstrip().split("\t")
        initialDict[line[manifestExternalIdIndex]] = [line[fileNameIndex], line[familyIdIndex], line[probandIndex]]

    biospecimenColumnNames = biospecimen.readline()
    biospecimenColumnNames = biospecimenColumnNames.rstrip().split("\t")
    biospecimenExternalIdIndex = biospecimenColumnNames.index("External Aliquot Id")
    sampleIdIndex = biospecimenColumnNames.index("Biospecimens Id")

    for line in biospecimen:
        line = line.rstrip().split("\t")
        initialDict[line[biospecimenExternalIdIndex]].append(line[sampleIdIndex])

    for key, value in initialDict.items():
        finalDict[value[0]] = [value[1], value[3], value[2]]

    output.write("fileName\tfamilyID\tsampleID\tproband\n")
    for key, value in finalDict.items():
        output.write("{}\t{}\t{}\t{}\n".format(key, value[0], value[1], value[2]))
    