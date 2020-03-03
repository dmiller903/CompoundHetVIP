import os
import time
import argparse
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description="Annotate VCF using snpEff")

parser.add_argument('input_vcf', help='Annotated VCF File')
parser.add_argument('output_database', help='Name of output database (name needs to end in .db)')
parser.add_argument('--fam_file', help="If you have family's in the VCF file, include a fam file")
parser.add_argument('--num_cores', help='Loading will go quicker if more cores are available.', default=2)

args = parser.parse_args()

#Create variables of each argument from argparse
inputFile = args.input_vcf
databaseName = args.output_database
famFile = args.fam_file
cores = args.num_cores

# Download annotation files and CADD files into the container
if len(os.listdir('/usr/local/share/gemini/gemini_data')) == 0:
    os.system("gemini update --dataonly --extra cadd_score")

# Load annotated file into a GEMINI database
os.system(f"gemini load -v {inputFile} \
-p {famFile} -t snpEff --cores {cores} {databaseName}")

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Annotation complete. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))