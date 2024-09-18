import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from sklearn.preprocessing import StandardScaler

def fetch_all_data(db_path):
    """Fetch all molecule data from the database."""
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT MoleculeData.name, MoleculeData.vector, 
                   GROUP_CONCAT(SpectrumData.wavelength) as wavelengths,
                   GROUP_CONCAT(SpectrumData.intensity) as intensities
            FROM MoleculeData
            JOIN SpectrumData ON MoleculeData.id = SpectrumData.molecule_id
            GROUP BY MoleculeData.id
        """)
        return cursor.fetchall()

def standardize_spectrum(wavelengths, intensities, standard_x):
    """Interpolate spectrum to standard x values and normalize intensities."""
    wavelengths = np.array([float(w) for w in wavelengths.split(',')])
    intensities = np.array([float(i) for i in intensities.split(',')])
    
    # Interpolate to standard x values
    f = interpolate.interp1d(wavelengths, intensities, kind='linear', fill_value='extrapolate')
    standardized_intensities = f(standard_x)
    
    # Normalize intensities
    standardized_intensities = (standardized_intensities - np.min(standardized_intensities)) / \
                               (np.max(standardized_intensities) - np.min(standardized_intensities))
    
    return standardized_intensities

def prepare_dataset(db_path, num_points=1000):
    data = fetch_all_data(db_path)
    
    if not data:
        print("No data fetched from the database. The database might be empty or there might be an issue with the query.")
        return None

    print(f"Fetched {len(data)} molecules from the database.")
    
    # Determine the global min and max wavelengths
    all_wavelengths = []
    for mol in data:
        try:
            all_wavelengths.extend([float(w) for w in mol[2].split(',')])
        except ValueError as e:
            print(f"Error processing wavelengths for molecule {mol[0]}: {e}")
            print(f"Problematic wavelength data: {mol[2]}")
            continue

    if not all_wavelengths:
        print("No valid wavelength data found.")
        return None

    min_wavelength, max_wavelength = min(all_wavelengths), max(all_wavelengths)
    print(f"Wavelength range: {min_wavelength} - {max_wavelength}")
    
    # Create standard x-axis
    standard_x = np.linspace(min_wavelength, max_wavelength, num_points)
    
    # Prepare data
    X = []  # Spectra
    y = []  # Functional group vectors
    names = []  # Molecule names
    
    for name, vector, wavelengths, intensities in data:
        try:
            standardized_intensities = standardize_spectrum(wavelengths, intensities, standard_x)
            X.append(standardized_intensities)
            y.append([int(v) for v in vector.split(',')])
            names.append(name)
        except Exception as e:
            print(f"Error processing molecule {name}: {e}")
            continue
    
    if not X:
        print("No valid spectra data processed.")
        return None

    X = np.array(X)
    y = np.array(y)
    
    # Standardize X
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    return X_scaled, y, names, standard_x

def plot_sample_spectra(X, names, standard_x, num_samples=5):
    plt.figure(figsize=(12, 6))
    for i in range(min(num_samples, len(X))):
        plt.plot(standard_x, X[i], label=names[i])
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('Standardized Intensity')
    plt.title('Sample Standardized Spectra')
    plt.legend()
    plt.tight_layout()
    plt.savefig('sample_standardized_spectra.png')
    plt.close()

def plot_functional_group_heatmap(y, names):
    functional_groups = [
        "Alkyl halide", "Alcohol", "Aldehyde", "Ketone", "Carboxylic acid",
        "Ester", "Ether", "Primary amine", "Secondary amine", "Tertiary amine",
        "Amide", "Nitro", "Nitrile", "Sulfide", "Sulfoxide", "Sulfone",
        "Phosphate", "Phenol", "Imine", "Alkene", "Alkyne", "Thiol",
        "Acyl chloride", "Anhydride", "Epoxide", "Lactam", "Aromatic", "Alicyclic"
    ]
    
    plt.figure(figsize=(15, 10))
    plt.imshow(y.T, aspect='auto', cmap='viridis')
    plt.colorbar(label='Presence of Functional Group')
    plt.xlabel('Molecule')
    plt.ylabel('Functional Group')
    plt.title('Functional Group Heatmap')
    plt.xticks(range(len(names)), names, rotation=90, ha='right')
    plt.yticks(range(len(functional_groups)), functional_groups)
    plt.tight_layout()
    plt.savefig('functional_group_heatmap.png')
    plt.close()

if __name__ == "__main__":
    db_path = 'labeled.db'
    result = prepare_dataset(db_path)
    
    if result is None:
        print("Failed to prepare dataset. Check the error messages above.")
    else:
        X, y, names, standard_x = result
        
        print(f"Dataset shape: X: {X.shape}, y: {y.shape}")
        print(f"Wavenumber range: {standard_x[0]:.2f} - {standard_x[-1]:.2f} cm⁻¹")
        
        # Plot some sample spectra
        plot_sample_spectra(X, names, standard_x)
        print("Sample standardized spectra plot saved as 'sample_standardized_spectra.png'")
        
        # Plot functional group heatmap
        plot_functional_group_heatmap(y, names)
        print("Functional group heatmap saved as 'functional_group_heatmap.png'")
