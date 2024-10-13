import os
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def parse_jdx(file_path):
    """Parse a JDX file and extract metadata, wavelengths, and intensities."""
    metadata = {}
    wavelengths = []
    intensities = []
    data_started = False

    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            
            if not data_started:
                # Extract metadata
                if line.startswith('##') and '=' in line:
                    key, value = line.split('=', 1)
                    metadata[key] = value
                elif line.startswith('##XYDATA='):
                    data_started = True
                continue  # Continue to next line until data starts

            # Data parsing
            if line.startswith('##END='):
                break  # End of data section
            if not line or line.startswith('##'):
                continue  # Skip empty lines or metadata within data section
            
            parts = line.split()
            if not parts:
                continue  # Skip empty lines

            # Attempt to parse the wavelength
            wavelength_str = parts[0]
            try:
                wavelength = float(wavelength_str)
            except ValueError:
                logging.warning(f"Line {line_number}: Unable to convert wavelength '{wavelength_str}' to float. Skipping line.")
                continue  # Skip lines with invalid wavelength

            # Parse intensity values
            for intensity_str in parts[1:]:
                try:
                    intensity = float(intensity_str)
                except ValueError:
                    logging.warning(f"Line {line_number}: Unable to convert intensity '{intensity_str}' to float. Skipping value.")
                    continue  # Skip invalid intensity values

                wavelengths.append(wavelength)
                intensities.append(intensity)

    # Convert lists to NumPy arrays
    wavelengths = np.array(wavelengths)
    intensities = np.array(intensities)

    return metadata, wavelengths, intensities

def plot_spectrum(wavenumbers, intensities, title="Infrared Spectrum"):
    plt.figure(figsize=(12, 8))
    if intensities.ndim == 1:
        plt.plot(wavenumbers, intensities, label="Spectrum")
    else:
        for i in range(intensities.shape[1]):
            plt.plot(wavenumbers, intensities[:, i], label=f"Spectrum {i+1}")
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('Transmittance')
    plt.title(title)
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    # Get all JDX files in the IR directory
    ir_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'IR')
    jdx_files = [f for f in os.listdir(ir_directory) if f.lower().endswith('.jdx')]
    if not jdx_files:
        raise FileNotFoundError("No .jdx files found in the IR directory.")

    logging.info(f"Found {len(jdx_files)} JDX files in '{ir_directory}' directory.")

    # Connect to SQLite database (or create it if it doesn't exist)
    db_path = 'spectrum_data.db'
    conn = sqlite3.connect(db_path)
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
            wavelength REAL,
            intensity REAL,
            FOREIGN KEY (molecule_id) REFERENCES molecules (id)
        )
    ''')

    conn.commit()
    logging.info(f"Connected to SQLite database at '{db_path}' and ensured tables exist.")

    # Process each JDX file
    for jdx_file in jdx_files:
        file_path = os.path.join(ir_directory, jdx_file)
        try:
            metadata, wavelengths, intensities = parse_jdx(file_path)
        except ValueError as ve:
            logging.error(f"Error parsing {jdx_file}: {ve}")
            continue  # Skip to the next file

        # Extract desired metadata
        title = metadata.get('##TITLE', os.path.splitext(jdx_file)[0])
        logging.info(f"Processing File: {title}")

        # Insert molecule name into molecules table
        molecule_name = os.path.splitext(jdx_file)[0]
        cursor.execute('INSERT INTO molecules (molecule_name) VALUES (?)', (molecule_name,))
        molecule_id = cursor.lastrowid

        # Insert data into the spectrum_data table
        data = [(molecule_id, wl, inten) for wl, inten in zip(wavelengths, intensities)]
        cursor.executemany('INSERT INTO spectrum_data (molecule_id, wavelength, intensity) VALUES (?, ?, ?)', data)
        conn.commit()
        logging.info(f"Inserted data for '{molecule_name}' into the database.")

    # Retrieve and print the contents of the tables
    cursor.execute('SELECT * FROM molecules LIMIT 5')
    molecule_rows = cursor.fetchall()

    print("\nMolecules Table (First 5 Entries):")
    print("ID | Molecule Name")
    for row in molecule_rows:
        print(row)

    cursor.execute('SELECT * FROM spectrum_data LIMIT 10')
    spectrum_rows = cursor.fetchall()

    print("\nSpectrum Data Table (First 10 Entries):")
    print("ID | Molecule ID | Wavelength | Intensity")
    for row in spectrum_rows:
        print(row)

    # Close the database connection
    conn.close()
    logging.info("Closed the database connection.")

    logging.info("All data has been processed and stored in the database.")

if __name__ == "__main__":
    main()