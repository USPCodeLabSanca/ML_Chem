import os
from multiprocessing import Pool

import sqlite3
import matplotlib.pyplot as plt
import numpy as np

def run_process(args):
    process, array = args
    os.system(f'python3 {process} {" ".join(map(str, array))}')


# Connect to the SQLite database
conn = sqlite3.connect('final.db')
cursor = conn.cursor()


# Get the total number of molecules
cursor.execute('SELECT COUNT(*) FROM Molecules')
n_molecules = cursor.fetchone()[0]


# Just get the number of molecules
print(f"Total number of molecules: {n_molecules}")





n = 12  # Number of workers
#n_molecules = 100


# Calculate the division of molecules into n groups
group_sizes = [n_molecules // n] * (n - 1)
group_sizes.append(n_molecules - sum(group_sizes))

print(f"Number of groups: {n}")
print("Molecule distribution:")
for i, size in enumerate(group_sizes, 1):
    print(f"n{i} = {size}")



irange = []
for i in range(len(group_sizes)):
    start = sum(group_sizes[:i]) + 1
    end = sum(group_sizes[:i+1])
    irange.append([start, end])


# List of tuples containing process name and array to pass
processes = []
i = 1
for slice in irange:
    slice.append(i)
    processes.append(('test1.py', slice))
    i += 1
print("Processes:")
print(processes)

# Create a pool of n worker processes
with Pool(processes=n) as pool:
    # Map the run_process function to the processes list
    pool.map(run_process, processes)

print("Joining DBs")

print("All processes completed.")

