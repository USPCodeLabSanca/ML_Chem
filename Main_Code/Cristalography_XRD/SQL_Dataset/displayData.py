import sqlite3
import os
import numpy as np

# Connect to the SQLite database in the parent directory
conn = sqlite3.connect('rruff.db')
cursor = conn.cursor()

# Check if the 'rruff_xdr' table exists
cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='rruff_xdr'")
if cursor.fetchone() is None:
    print("Error: The 'rruff_xdr' table does not exist in the database.")
else:
    # Display the table contents in the terminal
    print("\nDisplaying the contents of the 'rruff_xdr' table:")
    print("-----------------------------------------------")

    try:
        # Fetch all rows from the rruff_xdr table
        cursor.execute('SELECT id, label, y_values FROM rruff_xdr LIMIT 5')
        rows = cursor.fetchall()

        # Print the table headers
        print(f"{'ID':<5}{'Label':<50}{'Y Values (first 5 points)':<50}")
        print("-" * 105)

        # Print each row
        for row in rows:
            id, label, y_values_blob = row
            y_values = np.frombuffer(y_values_blob, dtype=np.float32)
            y_values_str = ', '.join(map(str, y_values[:5]))  # Convert first 5 points to string
            print(f"{id:<5}{label:<50}{y_values_str:<50}")

        print("\nNote: Y values are stored as binary blobs in the database.")
        print("They are converted to numpy arrays for display.")
    except sqlite3.OperationalError as e:
        print(f"Error: {e}")

import matplotlib.pyplot as plt

# Fetch the first item from the database
cursor.execute('SELECT label, y_values FROM rruff_xdr LIMIT 1')
first_item = cursor.fetchone()

if first_item:
    label, y_values_blob = first_item
    y_values = np.frombuffer(y_values_blob, dtype=np.float32)
    
    # Create x values (assuming they are indices)
    x_values = np.arange(len(y_values))
    
    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values)
    plt.title(f'First Item from rruff_xdr Table: {label}')
    plt.xlabel('Index')
    plt.ylabel('Y Values')
    plt.show()
else:
    print("No data found in the rruff_xdr table.")

# Close the connection
conn.close()

# Database structure:
# Table name: rruff_xdr
# Columns:
#   - id: INTEGER PRIMARY KEY AUTOINCREMENT
#   - label: TEXT
#   - y_values: BLOB (contains binary data of numpy array)
#
# The 'rruff_xdr' table stores spectral data for different samples.
# Each row represents a single spectrum with the following information:
#   - id: A unique identifier for each spectrum
#   - label: The name or identifier of the sample
#   - y_values: The spectral intensity values stored as a binary blob
#     (This is typically a numpy array converted to binary format)
#
# How y_values are obtained:
# The y_values are likely obtained from spectral measurements of samples.
# These values represent the intensity or absorbance at different wavelengths or frequencies.
# In the original data processing (probably in IngestData.py), the y_values were extracted from
