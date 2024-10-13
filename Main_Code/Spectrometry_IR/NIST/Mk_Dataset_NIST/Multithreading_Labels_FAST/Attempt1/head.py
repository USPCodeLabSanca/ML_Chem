import os
from multiprocessing import Pool

def run_process(args):
    process, array = args
    os.system(f'python3 {process} {" ".join(map(str, array))}')



# List of tuples containing process name and array to pass
processes = [
    ('p1.py', [1, 2, 3]),
    ('p1.py', [4, 5, 6]),
    ('p1.py', [7, 8, 9])
]

# Create a pool of 3 worker processes
with Pool(processes=3) as pool:
    # Map the run_process function to the processes list
    pool.map(run_process, processes)

print("All processes completed.")
