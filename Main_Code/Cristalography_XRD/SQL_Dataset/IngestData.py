# Finding the Data Directory
import os
import torch
import sqlite3
import datetime
import numpy as np
import matplotlib.pyplot as plt

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


def Create_DataPoint(file_path):
    x_values, y_values = get_data(file_path)

    if x_values is None or y_values is None:
        return None, None, None  # Return None for all values if get_data returns None

    if len(x_values) != 8501:
        return None, None, None  # Return None for all values if the file doesn't have 8501 data points

    # Extract the molecule name from the filename
    filename = os.path.basename(file_path)
    label = filename.split('__')[0]  # Get the part before the first '__'

    # Convert to tensor
    x_tensor = torch.from_numpy(np.array(x_values)).float()
    y_tensor = torch.from_numpy(np.array(y_values)).float()

    return x_tensor, y_tensor, label


def process_all_files_in_directory():
    x_tensor = None
    y_tensor = None
    all_labels = []
    file_order = []

    # Get the current working directory
    directory_path = os.getcwd()

    # Get all files with .txt extension in the current directory and sort them alphabetically by label
    txt_files = sorted([f for f in os.listdir(directory_path) if f.endswith('.txt')], key=lambda x: x.split('__')[0])

    for file_name in txt_files:
        file_path = os.path.join(directory_path, file_name)
        x, y, label = Create_DataPoint(file_path)

        if x is not None and y is not None and label is not None:
            # Add data to the main tensors
            if x_tensor is None:
                x_tensor = x.unsqueeze(0)
                y_tensor = y.unsqueeze(0)
            else:
                x_tensor = torch.cat((x_tensor, x.unsqueeze(0)), dim=0)
                y_tensor = torch.cat((y_tensor, y.unsqueeze(0)), dim=0)

            all_labels.append(label)
            file_order.append(file_name)

    return x_tensor, y_tensor, all_labels, file_order

# Usage
x_tensor, y_tensor, labels, file_order = process_all_files_in_directory()

# ---------// Creating Database!! ----------------------------------

parent_dir = os.path.dirname(os.getcwd())

# Create or connect to the SQLite database in the parent directory
conn = sqlite3.connect(os.path.join(parent_dir, 'rruff.db'))
cursor = conn.cursor()

# Create a table to store the data
cursor.execute('''
    CREATE TABLE IF NOT EXISTS rruff_xdr (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        file_order INTEGER,
        label TEXT,
        x_values BLOB,
        y_values BLOB,
        date_added DATE
    )
''')

# Insert data into the database
current_date = datetime.date.today().isoformat()
for i, (label, file_name) in enumerate(zip(labels, file_order)):
    x_values = x_tensor[i].numpy().tobytes()
    y_values = y_tensor[i].numpy().tobytes()
    cursor.execute('INSERT INTO rruff_xdr (file_order, label, x_values, y_values, date_added) VALUES (?, ?, ?, ?, ?)', 
                   (i, label, x_values, y_values, current_date))

# Commit changes and close the connection
conn.commit()
conn.close()

print(f"Data has been successfully stored in {os.path.join(parent_dir, 'rruff.db')}")

'''
Database structure explanation:
The 'rruff_xdr' table is designed to store spectral data for XRD (X-ray diffraction) measurements.
Each row in the table represents a single XRD spectrum with the following columns:
- id: A unique identifier for each spectrum, automatically incremented.
- file_order: An integer representing the order of the file in the original dataset (now in alphabetical order by label).
- label: The name or identifier of the sample, extracted from the original filename.
- x_values: The x-axis values of the spectrum (typically 2θ angles in XRD), stored as a binary blob.
- y_values: The y-axis values of the spectrum (typically intensity in XRD), stored as a binary blob.
- date_added: The date when the data was added to the database.

The x_values and y_values are stored as binary blobs to efficiently handle large arrays of data.
These blobs are created by converting numpy arrays to binary format using the tobytes() method.
When retrieving the data, these blobs need to be converted back to numpy arrays for analysis or plotting.

This structure allows for efficient storage and retrieval of XRD spectral data, enabling further
analysis and comparison of different mineral samples in the RRUFF database.
The data is now stored in alphabetical order by the label, maintaining the original file order.
'''
