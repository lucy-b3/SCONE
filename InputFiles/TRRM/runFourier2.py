import os
import subprocess
import re
import numpy as np

SCONE = "./../../Build/scone.out"
INPUT_FILE = "slab_3"
OUTPUT_FILE = "slab_3.m"
N = 100

os.makedirs("results", exist_ok=True)

# Updated regex pattern for 3D reshape
pattern = re.compile(
    r"Iterate_Error_Iterate_Error\s*=\s*reshape\(\[\s*(.*?)\s*\],\s*2,\s*2,\s*(\d+)\);",
    re.DOTALL
)

all_runs = []

for i in range(N):
    print(f"Running SCONE iteration {i+1}/{N} ...")
    subprocess.run([SCONE, INPUT_FILE,"--omp", "25"])
    
    run_out = f"results/output_{i}.m"
    os.rename(OUTPUT_FILE, run_out)
    
    with open(run_out, "r") as f:
        text = f.read()
    
    match = pattern.search(text)
    if match:
        data_str = match.group(1)
        values = np.fromstring(data_str.replace(",", " "), sep=" ")
        dim3 = int(match.group(2))
        #arr = values.reshape(2, 2, dim3)
        arr = values.reshape(2,dim3,2)
        #arr = values.reshape((2,2,dim3), order='F')
        all_runs.append(arr)
        #values = np.fromstring(data_str.replace(",", " "), sep=" ")
        #dim3 = int(match.group(2))
        #arr = values.reshape(2, 2, dim3, order='F')
        #real_part = arr[0, :, :]  # take only real component
        #all_runs.append(real_part)

    else:
        print(f"Warning: 'Iterate_Error_Iterate_Error' not found in {run_out}")

    print('ITERATION ' + str(i+1) + ' OF ' + str(N)) 
        
# Combine into one array across runs
# Result shape: (N, 2, 2, 100)
all_runs = np.stack(all_runs)
print("Combined array shape:", all_runs.shape)

# Save for later analysis
np.save("combined_iterate_error.npy", all_runs)
    


