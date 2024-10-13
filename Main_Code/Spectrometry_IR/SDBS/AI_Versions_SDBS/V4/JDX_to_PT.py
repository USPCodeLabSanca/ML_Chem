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

def plot_spectrum(wavenumbers, intensities, title="Infrared Spectrum"):
    plt.figure(figsize=(10, 6))
    plt.plot(wavenumbers, intensities, label=title)
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('Transmittance')
    plt.title(title)
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.legend()
    plt.show()

# Get all JDX files in the IR directory
ir_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'IR')
jdx_files = [f for f in os.listdir(ir_directory) if f.endswith('.jdx')]
if not jdx_files:
    raise FileNotFoundError("No .jdx files found in the IR directory.")

# Determine common wavelength range and number of points
all_wavelengths = []
for file in jdx_files:
    wl, _ = parse_jdx(os.path.join(ir_directory, file))
    all_wavelengths.append(wl)

# Find global min and max wavelengths
min_wavelength = max(wl.min() for wl in all_wavelengths)
max_wavelength = min(wl.max() for wl in all_wavelengths)

# Define common wavelength grid
standard_npoints = min(len(wl) for wl in all_wavelengths)
standard_wavelengths = np.linspace(min_wavelength, max_wavelength, standard_npoints)

# Connect to SQLite database (or create it if it doesn't exist)
conn = sqlite3.connect('spectrum_data.db')
cursor = conn.cursor()

# Create tables
cursor.execute('''
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_name TEXT
    )
''')

cursor.execute('''
    CREATE TABLE IF NOT EXISTS spectrum_data (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INTEGER,
        n INTEGER,
        wavelength REAL,
        intensity REAL,
        FOREIGN KEY (molecule_id) REFERENCES molecules (id)
    )
''')

# Process each JDX file
standard_wavenumbers = 10000 / standard_wavelengths
intensity_matrix = []

for jdx_file in jdx_files:
    wavelengths, intensities = parse_jdx(os.path.join(ir_directory, jdx_file))
    # Interpolate intensities to the standard wavelength grid
    interp_func = interpolate.interp1d(wavelengths, intensities, kind='linear', fill_value="extrapolate")
    standardized_intensities = interp_func(standard_wavelengths)
    intensity_matrix.append(standardized_intensities)
    
    # Insert molecule name into molecules table
    molecule_name = os.path.splitext(jdx_file)[0]
    cursor.execute('INSERT INTO molecules (molecule_name) VALUES (?)', (molecule_name,))
    molecule_id = cursor.lastrowid
    
    # Insert data into the spectrum_data table
    data = [(molecule_id, n, wl, inten) for n, (wl, inten) in enumerate(zip(standard_wavelengths, standardized_intensities), 1)]
    cursor.executemany('INSERT INTO spectrum_data (molecule_id, n, wavelength, intensity) VALUES (?, ?, ?, ?)', data)
    conn.commit()

print(f"Processed {len(jdx_files)} JDX files with {standard_npoints} standardized points each.")

# Retrieve and print the contents of the tables
cursor.execute('SELECT * FROM molecules LIMIT 5')
molecule_rows = cursor.fetchall()

print("Molecules Table:")
print("ID | Molecule Name")
for row in molecule_rows:
    print(row)

cursor.execute('SELECT * FROM spectrum_data LIMIT 10')
spectrum_rows = cursor.fetchall()

print("\nSpectrum Data Table:")
print("ID | Molecule ID | n | Wavelength | Intensity")
for row in spectrum_rows:
    print(row)

# Close the database connection
conn.close()

# Plotting all standardized spectra
intensity_matrix = np.array(intensity_matrix)
plot_spectrum(standard_wavenumbers, intensity_matrix.T, title="Standardized IR Spectra")
