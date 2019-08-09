import gzip
import glob
import re
import os
import time

startTime = time.time()

# Keeps metadata and removes all non-variant sites
for file in glob.glob('*.g.vcf.gz'):
    fileName = re.findall(r'(.+)\.g\.vcf\.gz', file)[0]
    outputName = '{}_variants_only.vcf'.format(fileName)
    with gzip.open(file, 'rt') as gVCF, open(outputName, 'wt') as VCF:
        for line in gVCF:
            if line.startswith('#'):
                VCF.write(line)
            else:
                break
    os.system('zgrep -Ev "#|END" {} >> {}'.format(file, outputName))

# VT tools will split, left-align, and trim variants in order to reduce false negatives
for file in glob.glob('*only.vcf'):
    fileName = re.findall(r'(.+_only)\.vcf', file)[0]
    outputName = '{}_variants_only.vcf'.format(fileName)
    os.system('vt decompose -s /proj/{} | vt normalize -n -r /references/Homo_sapiens_assembly38.fasta - > /proj/{}'.format(file, outputName))

# Annotate vcf's with snpEff
for file in glob.glob('*_vt.vcf'):
    fileName = re.findall(r'(.+)_variants.+', file)[0]
    outputName = '{}_annotated.vcf'.format(fileName)
    os.system('java -Xmx40g -jar /opt/conda/share/snpeff-4.3.1t-2/snpEff.jar -v GRCh37.75 \
        /proj/{} > /proj/{}'.format(file, outputName))

# Load vcf's into GEMINI to create GEMINI databases
for file in glob.glob('*_annotated.vcf'):
    fileName = re.findall('(.+)_annotated.vcf', file)[0]
    outputName = '{}.db'.format(fileName)
    os.system('gemini load -v /proj/{} -t snpEff --cores 40 /proj/{}'.format(file, outputName))

print('Done! Processed in {} minutes'.format(round((time.time()-startTime)/60)))

    
