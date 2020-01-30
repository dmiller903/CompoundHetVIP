import re
import os
from sys import argv

#Get disease name based on path
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

diseaseName = chromosome = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]

# Use VT to split, trim and left align the phased samples.
os.system(f"/root/miniconda2/bin/vt decompose -s {pathToFiles}/{diseaseName}_phased_samples.vcf.gz \
| /root/miniconda2/bin/vt normalize -n -r /references/human_g1k_v37.fasta - > \
{pathToFiles}/{diseaseName}_phased_samples_vt.vcf")