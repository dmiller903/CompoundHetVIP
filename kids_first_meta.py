#This function takes 4 arguments: manifest file, biospeciment file, clinical file, and the name/location of output file.

#Import necessary packages
import time
import argparse
import pandas as pd

#Save time that script started
startTime = time.time()

#Input/output files
parser = argparse.ArgumentParser(
    description="Combine a manifest, biospecimen, and clinical file into an output TSV"
)
parser.add_argument(
    "manifestFile",
    help='Path to a manifest TSV file with columns "File Name", "Family Id", "Aliquot External ID", and "Proband"'
)
parser.add_argument(
    "biospecimenFile",
    help='Path to a biospecimen TSV file with columns "External Aliquot Id" and "Biospecimens Id"'
)
parser.add_argument(
    "clinicalFile",
    help='Path to a clinical TSV file with columns "External Id" and "Gender"'
)
parser.add_argument(
    "outputFile",
    help='Path to the output file'
)
args = parser.parse_args()

# Load TSV files
biospecimen_df = pd.read_csv(
    args.biospecimenFile,
    sep="\t",
    usecols=["External Aliquot Id", "Biospecimens Id"]
).rename(columns={
    "External Aliquot Id": "external_id",
    "Biospecimens Id": "sample_id",
})
clinical_df = pd.read_csv(
    args.clinicalFile,
    sep="\t",
    usecols=["External Id", "Gender"]
).rename(columns={
    "External Id": "external_id",
    "Gender": "sex",
})
manifest_df = pd.read_csv(
    args.manifestFile,
    sep="\t",
    usecols=["Aliquot External ID", "File Name", "Family Id", "Proband"]
).rename(columns={
    "Aliquot External ID": "external_id",
    "File Name": "file_name",
    "Family Id": "family_id",
    "Proband": "proband",
})

# Join the dataframes
df = manifest_df.merge(biospecimen_df)\
    .merge(clinical_df)\
    .drop(columns=["external_id"])

# Keep families with all 3 members
keep = df.groupby(["family_id"]).size() == 3
df = df[df["family_id"].isin(keep[keep].index)] \
    .sort_values(["family_id"])

# Reorder columns
df = df[["file_name", "family_id", "sample_id", "proband", "sex"]]

# Convert sex to numeric values
def sex_to_number(sex):
    if sex == "Female":
        return 2
    return 1
df["sex"] = df["sex"].apply(sex_to_number)

# Save output file
df.to_csv(args.outputFile, sep="\t", index=False)

#output message
timeElapsedSeconds = round((time.time()-startTime), 2)
outputMessage = "Combined manifest, biospecimen, and clinical files. Time elapsed: {} seconds".format(timeElapsedSeconds)
lenMessage = len(outputMessage)
char = '\n' + ('*' * lenMessage) + '\n'
print('{}{}{}'.format(char, outputMessage, char))
