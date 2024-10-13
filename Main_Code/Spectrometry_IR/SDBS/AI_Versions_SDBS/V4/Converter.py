import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import torch
import sqlite3

def parse_jdx(file_path):
    """Parse a JDX file and extract wavelengths and intensities."""
    wavelengths = []
    intensities = []
    deltax = None
    firstx = None
    lastx = None
    npoints = None

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('##DELTAX='):
                deltax = float(line.split('=')[1])
            elif line.startswith('##FIRSTX='):
                firstx = float(line.split('=')[1])
            elif line.startswith('##LASTX='):
                lastx = float(line.split('=')[1])
            elif line.startswith('##NPOINTS='):
                npoints = int(line.split('=')[1])
            elif line.startswith('##XYDATA='):
                break

        if None in (deltax, firstx, lastx, npoints):
            raise ValueError("Missing required metadata in JDX file")

        data_started = False
        current_wavelength = firstx
        for line in file:
            if line.startswith('##END='):
                break
            if not data_started:
                data_started = True
                continue
            parts = line.strip().split()
            if not parts:
                continue
            wavelength = float(parts[0])
            for intensity in parts[1:]:
                wavelengths.append(wavelength)
                intensities.append(float(intensity))
                wavelength += deltax

    return np.array(wavelengths), np.array(intensities)


def plot_spectrum(wavelengths, intensities, title="Infrared Spectrum"):
    # Convert wavelengths in micrometers to wavenumbers in cm⁻¹
    wavenumbers = 10000 / wavelengths
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wavenumbers, intensities, color='blue')
    ax.set_xlabel('Wavenumber (cm⁻¹)')
    ax.set_ylabel('Transmittance')
    ax.set_title(title)
    ax.invert_xaxis()
    ax.grid(True)
    plt.show()

# Get the first JDX file in the current directory
jdx_files = [f for f in os.listdir() if f.endswith('.jdx')]
if jdx_files:
    jdx_file = jdx_files[0]
else:
    raise FileNotFoundError("No .jdx files found in the current directory.")

# Parse the JDX file
wavelengths, intensities = parse_jdx(jdx_file)

# Print the number of points
print(f"Number of data points: {len(wavelengths)}")

# Connect to SQLite database (or create it if it doesn't exist)
conn = sqlite3.connect('spectrum_data.db')
cursor = conn.cursor()

# Create table
cursor.execute('''
    CREATE TABLE IF NOT EXISTS spectrum (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        wavelength REAL,
        intensity REAL
    )
''')

# Insert data into the table
data = list(zip(wavelengths, intensities))
cursor.executemany('INSERT INTO spectrum (wavelength, intensity) VALUES (?, ?)', data)
conn.commit()

# Retrieve and print the contents of the table
cursor.execute('SELECT * FROM spectrum')
rows = cursor.fetchall()

print("Spectrum Data Table:")
print("ID | Wavelength | Intensity")
for row in rows[:10]:  # Print only the first 10 entries for an overview
    print(row)

# Close the database connection
conn.close()

# plot_spectrum(wavelengths, intensities, title=f"IR Spectrum of {jdx_file[:-4]}")
