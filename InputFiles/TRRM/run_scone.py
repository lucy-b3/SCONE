#Runs SCONE via python

import subprocess, os
import numpy as np
import json


SCONE = "./../../build/scone.out"
#input = "thermalPin"
#output = "thermalPinOut"
input = "kaistPin"
output = "kaistPinOut"
outputFormat = ".json"

subprocess.call([SCONE,input]) 

fileName = output+outputFormat

"""

with open(fileName, 'r') as file:
    #print(type(file))
    #print(file.readlines())
    #for line in file.readlines():
        #print(line, end='')
    lines = file.readlines()
    
with open(fileName,'a') as file:
    for line in reversed(lines):
        file.write(line)
    
with open(fileName, 'r') as file:
    for line in file.readlines():
        print(line, end='')


with open(fileName, 'r') as file:
    for line in file:
        line = line.rstrip()
        if "flux" in line:
            print(line)
            
"""         
   
with open(fileName) as file:
    d = json.load(file) #imports into dictionary
    #print(d)
    
flux = np.array(d["flux1G"]["flux1G"])
#print(np.shape(power))
flux1D = np.squeeze(flux)
print("Flux values = ", flux1D[:,0])
print("Flux standard deviation = ", flux1D[:,1])

power = np.array(d["power"]["power"])
#print(np.shape(power))
power1D = np.squeeze(power)
print("Power values = ", power1D[:,0])
print("Power standard deviation = ", power1D[:,1])


temperature = np.linspace(700,800,20)

temp_list = " ".join(str(t) for t in temperature)
new_line = f"             fuel ( {temp_list} );"

updated_lines = []
with open(input, "r") as f:
    for line in f:
        stripped = line.strip()
        if stripped.startswith("fuel ("):
            updated_lines.append(new_line + "\n")
        else:
           updated_lines.append(line)

with open(input, "w") as f:
    f.writelines(updated_lines)


#subprocess.call([SCONE,input]) 


# Run WIMS-ARTHUR to generate output file
run_cmd = SCONE + " " + input

"""
with subprocess.Popen(run_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT) as p:
    # Capture all WIMS-ARTHUR stdout to logger, possible as wait called next (i.e. start/stop process)
    log_subprocess_output(logger,p.stdout)
        
# for WIMS to finish within its Python subprocess
p.wait()
"""
process = subprocess.Popen(run_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
stdout, stderr = process.communicate()

print(stdout)
print(' ')
print(stderr)
