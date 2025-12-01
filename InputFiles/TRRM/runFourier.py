import os
import subprocess
import json
import numpy as np

# === CONFIGURATION ===
SCONE = "./../../Build/scone.out"
INPUT_FILE = "slab_3"
OUTPUT_FILE = "slab_3.json"   # JSON output now
N = 25                        # Number of runs
os.makedirs("results", exist_ok=True)

# === STORAGE FOR ALL RUNS ===
all_runs = []

# === LOOP OVER RUNS ===
for i in range(N):
    print(f"Running SCONE iteration {i+1}/{N} ...")

    # Run SCONE (Fortran executable)
    subprocess.run([SCONE, INPUT_FILE, "--omp", "25"])

    # Move output file to a results folder
    run_out = f"results/output_{i}.json"
    os.rename(OUTPUT_FILE, run_out)

    # Load the JSON file
    with open(run_out, "r") as f:
        data = json.load(f)

    # Extract the "Iterate_Error" block
    iterate_error = data.get("Iterate_Error", {}).get("Iterate_Error", None)

    if iterate_error is None:
        print(f"Warning: 'Iterate_Error' not found in {run_out}")
        continue

    # Convert to NumPy array
    arr = np.array(iterate_error, dtype=float)
    arr = arr.transpose((1, 0, 2))   # (2, 500, 2)
    arr = arr[::-1, :, :]            # reverse iteration axis (swap j=1 and j=2)
    print("arr.shape:", arr.shape)
    print("arr[0, :, 0] first cell two iterations:", arr[0, :, 0])
    print("arr[1, :, 0] second cell two iterations:", arr[1, :, 0])
    print("arr.transpose((1,0,2))[0,0,:] after transpose:", arr.transpose((1,0,2))[0,0,:])
    # Sanity check: arr should have shape (nCells-1, 2, 2)
    # Your Fortran writes the array as [i, j, 2]
    # You can reorder if needed:
    #arr = np.transpose(arr, (1, 0, 2))  # Now shape = (2, nCells-1, 2)
    # Determine shape automatically
    #if arr.shape[0] == 2 and arr.shape[1] != 2:
        # Already in correct shape (iteration, cells, value_type)
     #   pass
    #elif arr.shape[1] == 2 and arr.shape[0] != 2:
        # Needs transposing (was cells × iteration × value_type)
     #   arr = np.transpose(arr, (1, 0, 2))
    #else:
     #   print(f"⚠️ Unexpected shape {arr.shape} for Iterate_Error")
    print(arr.shape)
    arr = arr.transpose((1, 0, 2))  # now shape (2, 500, 2)
    all_runs.append(arr)

    print(f"ITERATION {i+1} OF {N} complete, shape {arr.shape}")

# === COMBINE ALL RUNS ===
# Stack along a new axis: (N_runs, 2, nCells-1, 2)
all_runs = np.stack(all_runs)
print("Combined array shape:", all_runs.shape)

# === SAVE ===
np.save("combined_iterate_error.npy", all_runs)
