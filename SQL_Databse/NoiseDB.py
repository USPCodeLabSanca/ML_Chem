# Finding the Data Directory
import os
import torch
import sqlite3
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# -------------------------------------------

print(os.listdir())
os.chdir('XY_Processed')
print(os.listdir())

print("------------")

dirs = os.listdir()

# -------------------------------------------

# Get Data Function
def get_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    metadata = []
    for line in lines:
        if line.startswith('##'):
            metadata.append(line.strip())
        else:
            try:
                x, y = map(float, line.strip().split(','))
                data.append((x, y))
            except ValueError:
                pass

    if not data:
        print("Error: No valid data points found in the file.")
        return

    x_values, y_values = zip(*data)
    return x_values, y_values

# -------------------------------------------

# Creating more Examples for each Sample by adding noise and X-shift
def create_noise(x_values, y_values, n=5, intensity_noise=0.04, peak_shift=1, peak_threshold=0.5):
    x_values = np.array(x_values)
    y_values = np.array(y_values)

    # Identify peaks
    peaks, _ = find_peaks(y_values, height=peak_threshold * np.max(y_values), distance=1)

    noisy_data = []

    for _ in range(n):
        # Create peak shifts
        shifts = np.random.uniform(-peak_shift, peak_shift, size=len(peaks))
        shift_array = np.zeros_like(x_values)
        shift_array[peaks] = shifts

        # Smooth out shifts for neighboring points
        kernel = np.array([0.5, 1, 0.5])
        smoothed_shifts = np.convolve(shift_array, kernel, mode='same') / np.sum(kernel)

        noisy_x = x_values + smoothed_shifts

        # Add intensity noise
        noise = np.random.uniform(-intensity_noise * np.max(y_values),
                                  intensity_noise * np.max(y_values),
                                  size=len(y_values))
        noisy_y = np.maximum(0, y_values + noise)

        # Sort the data points based on x values to maintain order
        sorted_indices = np.argsort(noisy_x)
        noisy_x = noisy_x[sorted_indices]
        noisy_y = noisy_y[sorted_indices]

        noisy_data.append((noisy_x, noisy_y))

    return noisy_data

def Create_DataPoint(file_path, n=5):
    x_values, y_values = get_data(file_path)

    if x_values is None or y_values is None:
        return None, None, None, None, None  # Return None for all values if get_data returns None

    if len(y_values) != 8501:
        return None, None, None, None, None  # Return None for all values if the file doesn't have 8501 data points

    # Extract the molecule name from the filename
    filename = os.path.basename(file_path)
    label = filename.split('__')[0]  # Get the part before the first '__'

    # Generate noisy data
    noisy_data = create_noise(x_values, y_values, n)

    # Convert to tensor
    x_tensor = torch.from_numpy(np.array(x_values)).float()
    y_tensor = torch.from_numpy(np.array(y_values)).float()

    # Convert noisy data to tensors
    noisy_x_tensors = [torch.from_numpy(np.array(nx)).float() for nx, _ in noisy_data]
    noisy_y_tensors = [torch.from_numpy(np.array(ny)).float() for _, ny in noisy_data]

    return x_tensor, y_tensor, noisy_x_tensors, noisy_y_tensors, label

def process_all_files_in_directory(n=5):
    x_tensor = None
    y_tensor = None
    noisy_x_tensors = []
    noisy_y_tensors = []
    all_labels = []
    file_order = []

    # Get the current working directory
    directory_path = os.getcwd()

    # Get all files with .txt extension in the current directory and sort them alphabetically by label
    txt_files = sorted([f for f in os.listdir(directory_path) if f.endswith('.txt')], key=lambda x: x.split('__')[0])

    for file_name in txt_files:
        file_path = os.path.join(directory_path, file_name)
        x, y, nx, ny, label = Create_DataPoint(file_path, n)

        if x is not None and y is not None and label is not None:
            # Add data to the main tensors
            if x_tensor is None:
                x_tensor = x.unsqueeze(0)
                y_tensor = y.unsqueeze(0)
            else:
                x_tensor = torch.cat((x_tensor, x.unsqueeze(0)), dim=0)
                y_tensor = torch.cat((y_tensor, y.unsqueeze(0)), dim=0)

            noisy_x_tensors.extend(nx)
            noisy_y_tensors.extend(ny)
            all_labels.extend([label] * (n + 1))
            file_order.extend([file_name] * (n + 1))

    return x_tensor, y_tensor, noisy_x_tensors, noisy_y_tensors, all_labels, file_order

# Usage
x_tensor, y_tensor, noisy_x_tensors, noisy_y_tensors, labels, file_order = process_all_files_in_directory(n=5)

# ---------// Creating Database!! ----------------------------------

parent_dir = os.path.dirname(os.getcwd())

# Create or connect to the SQLite database in the parent directory
conn = sqlite3.connect(os.path.join(parent_dir, 'rruff_noised.db'))
cursor = conn.cursor()

# Create a table to store the data
cursor.execute('''
    CREATE TABLE IF NOT EXISTS rruff_xdr (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        file_order INTEGER,
        label TEXT,
        x_values BLOB,
        y_values BLOB,
        is_noisy BOOLEAN,
        date_added DATE
    )
''')

# Insert data into the database
current_date = datetime.date.today().isoformat()
for i, (label, file_name) in enumerate(zip(labels, file_order)):
    if i % 6 == 0:  # Original data
        x_values = x_tensor[i//6].numpy().tobytes()
        y_values = y_tensor[i//6].numpy().tobytes()
        is_noisy = False
    else:  # Noisy data
        x_values = noisy_x_tensors[i-1].numpy().tobytes()
        y_values = noisy_y_tensors[i-1].numpy().tobytes()
        is_noisy = True
    
    cursor.execute('INSERT INTO rruff_xdr (file_order, label, x_values, y_values, is_noisy, date_added) VALUES (?, ?, ?, ?, ?, ?)', 
                   (i, label, x_values, y_values, is_noisy, current_date))

# Commit changes and close the connection
conn.commit()
conn.close()

print(f"Data has been successfully stored in {os.path.join(parent_dir, 'rruff_noised.db')}")

'''
Database structure explanation:
The 'rruff_xdr' table is designed to store spectral data for XRD (X-ray diffraction) measurements, including both original and noisy versions.
Each row in the table represents a single XRD spectrum with the following columns:
- id: A unique identifier for each spectrum, automatically incremented.
- file_order: An integer representing the order of the file in the dataset.
- label: The name or identifier of the sample, extracted from the original filename.
- x_values: The x-axis values of the spectrum (typically 2θ angles in XRD), stored as a binary blob.
- y_values: The y-axis values of the spectrum (typically intensity in XRD), stored as a binary blob.
- is_noisy: A boolean flag indicating whether the spectrum is an original (False) or a noisy version (True).
- date_added: The date when the data was added to the database.

The x_values and y_values are stored as binary blobs to efficiently handle large arrays of data.
These blobs are created by converting numpy arrays to binary format using the tobytes() method.
When retrieving the data, these blobs need to be converted back to numpy arrays for analysis or plotting.

This structure allows for efficient storage and retrieval of both original and noisy XRD spectral data, enabling further
analysis and comparison of different mineral samples in the RRUFF database, as well as providing augmented data for machine learning purposes.
The data is stored in the order of original sample followed by its noisy versions, maintaining the original file order.
'''
