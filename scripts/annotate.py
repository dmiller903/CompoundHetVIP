import re
import os
from sys import argv

#Get disease name based on path
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

diseaseName = chromosome = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]

# Annotate the vt trimmed file
os.system(f"java -Xmx40g -jar /snpEff/snpEff.jar GRCh37.75 -v \
{pathToFiles}/{diseaseName}_phased_samples_vt.vcf > \
{pathToFiles}/{diseaseName}_phased_samples_annotated.vcf")