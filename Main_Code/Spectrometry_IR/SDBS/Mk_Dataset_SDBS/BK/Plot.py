import sqlite3
import matplotlib.pyplot as plt
import numpy as np

def fetch_molecule_data(db_path, limit=5):
    """Fetch molecule data from the database."""
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT MoleculeData.name, MoleculeData.vector, 
                   GROUP_CONCAT(SpectrumData.wavelength) as wavelengths, 
                   GROUP_CONCAT(SpectrumData.intensity) as intensities
            FROM MoleculeData
            JOIN SpectrumData ON MoleculeData.id = SpectrumData.molecule_id
            GROUP BY MoleculeData.id
            LIMIT ?
        """, (limit,))
        return cursor.fetchall()

def plot_spectra(molecules):
    """Plot IR spectra for the given molecules."""
    plt.figure(figsize=(12, 8))
    
    for name, vector, wavelengths, intensities in molecules:
        wavelengths = [float(w) for w in wavelengths.split(',')]
        intensities = [float(i) for i in intensities.split(',')]
        
        plt.plot(wavelengths, intensities, label=name)
    
    plt.xlabel('Wavelength (cm⁻¹)')
    plt.ylabel('Intensity')
    plt.title('IR Spectra')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_functional_group_heatmap(molecules):
    """Plot a heatmap of functional groups for the given molecules."""
    names = [mol[0] for mol in molecules]
    vectors = [list(map(int, mol[1].split(','))) for mol in molecules]
    
    functional_groups = [
        "Alkyl halide", "Alcohol", "Aldehyde", "Ketone", "Carboxylic acid",
        "Ester", "Ether", "Primary amine", "Secondary amine", "Tertiary amine",
        "Amide", "Nitro", "Nitrile", "Sulfide", "Sulfoxide", "Sulfone",
        "Phosphate", "Phenol", "Imine", "Alkene", "Alkyne", "Thiol",
        "Acyl chloride", "Anhydride", "Epoxide", "Lactam", "Aromatic", "Alicyclic"
    ]
    
    plt.figure(figsize=(15, 10))
    plt.imshow(np.array(vectors).T, aspect='auto', cmap='viridis')
    plt.colorbar(label='Presence of Functional Group')
    plt.xlabel('Molecule')
    plt.ylabel('Functional Group')
    plt.title('Functional Group Heatmap')
    plt.xticks(range(len(names)), names, rotation=45, ha='right')
    plt.yticks(range(len(functional_groups)), functional_groups)
    plt.tight_layout()
    plt.show()

def get_wavelength_range(db_path):
    """Get the range of all wavelengths in the database."""
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT MIN(wavelength), MAX(wavelength)
            FROM SpectrumData
        """)
        min_wavelength, max_wavelength = cursor.fetchone()
    return min_wavelength, max_wavelength

if __name__ == "__main__":
    db_path = 'labeled.db'
    molecules = fetch_molecule_data(db_path, limit=5)
    
    if molecules:
        plot_spectra(molecules)
        plot_functional_group_heatmap(molecules)
        
        min_wavelength, max_wavelength = get_wavelength_range(db_path)
        print(f"Wavelength range in the database: {min_wavelength:.2f} - {max_wavelength:.2f} cm⁻¹")
    else:
        print("No data found in the database.")
