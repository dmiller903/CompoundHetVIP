import re
import os
import glob

#Download original files
os.system("wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz && tar xzf 1000GP_Phase3.tgz")

#Update original genetic map files format to Eagle and Beagle formats
for file in glob.glob("/references/1000GP_Phase3/genetic_map_chr*_combined_b37.txt"):
    fileName = re.findall(r"([\w_/]+)\.txt", file)[0]
    eagleOutput = f"{fileName}_eagle.txt"
    beagleOutput = f"{fileName}_beagle.txt"
    with open(file) as inputFile, open(eagleOutput, 'w') as eagleOutput, open(beagleOutput, 'w') as beagleOutput:
        header = inputFile.readline()
        header = "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
        eagleOutput.write(header)
        for line in inputFile:
            #Eagle format
            eagleLine = "22 " + line
            eagleOutput.write(eagleLine)
            #Beagle format
            beagleLineList = line.rstrip().split(" ")
            beagleLine = f"22  {beagleLineList[1]} {beagleLineList[2]} {beagleLineList[0]}\n"
            beagleOutput.write(beagleLine)