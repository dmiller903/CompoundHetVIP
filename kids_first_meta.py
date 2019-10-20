#This function takes 4 arguments: manifest file, biospeciment file, clinical file, and the name/location of output file.

#Import necessary packages
import time
from sys import argv

#Save time that script started
startTime = time.time()

#Input/output files
manifestFile = argv[1]
biospecimenFile = argv[2]
clinicalFile = argv[3]
outputFile = argv[4]

#Dictionary to store needed information in
outputDict = {}

#Obtain information from manifest file and add to initial dictionary
with open(manifestFile) as manifest:
    manifestColumnNames = manifest.readline()
    manifestColumnNames = manifestColumnNames.rstrip().split("\t")
    #Information needed from this file are the "File Name", "Family Id", "Aliquot External ID", and "Proband" (Yes, or No).
    fileNameIndex = manifestColumnNames.index("File Name")
    familyIdIndex = manifestColumnNames.index("Family Id")
    manifestExternalIdIndex = manifestColumnNames.index("Aliquot External ID")
    probandIndex = manifestColumnNames.index("Proband")
    
    #Add to the initial dictionary where the key is the "Aliquot External ID" as this is common across all input files.
    #The value is a list where value[0] the "File Name", value[1] is "Family Id", and value[2] is "Proband"
    #as "Yes" or "No"
    for sample in manifest:
        sample = sample.rstrip().split("\t")
        outputDict[sample[manifestExternalIdIndex]] = [sample[fileNameIndex], sample[familyIdIndex], sample[probandIndex]]

#Obtain information from the biospecimen file and add to initial dictionary
with open(biospecimenFile) as biospecimen:
    biospecimenColumnNames = biospecimen.readline()
    biospecimenColumnNames = biospecimenColumnNames.rstrip().split("\t")
    #Information needed from this file are the "External Aliquot Id", and "Biospecimens Id"
    biospecimenExternalIdIndex = biospecimenColumnNames.index("External Aliquot Id")
    sampleIdIndex = biospecimenColumnNames.index("Biospecimens Id")

    #Use the patient "External Aliquot Id" as the key to append the "Biospecimens Id" to the value list.
    for line in biospecimen:
        line = line.rstrip().split("\t")
        outputDict[line[biospecimenExternalIdIndex]].append(line[sampleIdIndex])

#Obtain gender information from the clinical file and add to initial dictionary
with open(clinicalFile) as clinical:
    clinicalColumnNames = clinical.readline()
    clinicalColumnNames = clinicalColumnNames.rstrip().split("\t")
    #Information needed from this file are the "External Id", and "Gender"
    clinicalExternalIdIndex = clinicalColumnNames.index("External Id")
    genderIndex = clinicalColumnNames.index("Gender")

    #Use the patient "External Id" as the key to append the "Biospecimens Id" to the value list.
    for line in clinical:
        line = line.rstrip().split("\t")
        outputDict[line[clinicalExternalIdIndex]].append(line[genderIndex])

#Create list of samples where only trios are included
familyDict = {}
for key, value in outputDict.items():
    if value[1] not in familyDict:
        familyDict[value[1]] = [value]
    else:
        familyDict[value[1]].append(value)

trioList = []
for key, value in familyDict.items():
    if len(value) == 3:
        trioList.append(value[0])
        trioList.append(value[1])
        trioList.append(value[2])

#Output information to outputFile
with open(outputFile, 'w') as output:
    output.write("file_name\tfamily_id\tsample_id\tproband\tsex\n")
    for item in trioList:
        if item[-1] == "Female":
            output.write("{}\t{}\t{}\t{}\t2\n".format(item[0], item[1], item[3], item[2], item[4]))
        else:
            output.write("{}\t{}\t{}\t{}\t1\n".format(item[0], item[1], item[3], item[2], item[4]))

#output message
timeElapsedSeconds = round((time.time()-startTime), 2)
outputMessage = "Combined manifest, biospecimen, and clinical files. Time elapsed: {} seconds".format(timeElapsedSeconds)
lenMessage = len(outputMessage)
char = '\n' + ('*' * lenMessage) + '\n'
print('{}{}{}'.format(char, outputMessage, char))