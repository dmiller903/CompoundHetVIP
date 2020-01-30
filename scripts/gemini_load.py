import re
import os
from sys import argv

#Get disease name based on path
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

diseaseName = chromosome = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]

# Load annotated file into a GEMINI database
os.system(f"gemini load -v {pathToFiles}/{diseaseName}_phased_samples_annotated.vcf \
-p {pathToFiles}/{diseaseName}.fam -t snpEff --cores 42 {pathToFiles}/{diseaseName}_phased_samples_annotated_cadd.db")