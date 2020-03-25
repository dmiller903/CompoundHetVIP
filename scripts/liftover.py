# Import necessary modules
import os
import time
import argparse
import re
import glob

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Some phasing programs and other programs such as GEMINI, require that VCF \
 or gVCF files be aligned to GRCh37. Therefore, this step takes an input VCF/gVCF that is aligned to GRCh38 converts it \
to GRCh37 using GATKs "CombineGVCFs" tool')

parser.add_argument('input_vcf', help='Input VCF file')
parser.add_argument('output_vcf', help='Path and name of output VCF file. Include ".gz" at the end of the file \
name')

args = parser.parse_args()

# Create variables of each argument from argparse
inputFile = args.input_vcf
outputFile = args.output_vcf
fileWithoutSuffix = re.findall(r'([\w\-_/]+)\.', outputFile)[0]

# Download necessary reference files if needed
if not os.path.exists("/references/hg38ToHg19.over.chain"):
    os.system("wget --no-check-certificate \
        https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d9ddd0ba7bc73000ce87e38/?zip= -O /tmp/references.zip \
        && unzip /tmp/references.zip -d /tmp/references \
        && rm /tmp/references.zip \
        && /root/miniconda2/bin/bgzip -d /tmp/references/human_g1k_v37_modified.fasta.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/hg38ToHg19.over.chain.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/Homo_sapiens_assembly38.fasta.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/Homo_sapiens_assembly38.dict.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/Homo_sapiens_assembly38.fasta.fai.gz")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}"):
            os.system(f"mv {file} /references/")
        elif fileName == "readme":
            os.system(f"mv {file} /references/")
    os.system("chmod 777 /references/*")

# Liftover file(s)
os.system(f"java -Xmx12G -jar /root/miniconda2/share/picard-2.21.1-0/picard.jar LiftoverVcf I={inputFile} O={outputFile} \
CHAIN=/references/hg38ToHg19.over.chain R=/references/human_g1k_v37_modified.fasta REJECT={fileWithoutSuffix}_rejected_variants.vcf")

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')