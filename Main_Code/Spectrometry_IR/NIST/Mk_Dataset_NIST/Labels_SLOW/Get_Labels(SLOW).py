import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem
import requests

def plot_spectrum(wavelengths, intensities, title="Infrared Spectrum"):
    wavelengths = np.array(wavelengths)
    wavenumbers = wavelengths
    intensities = np.array(intensities)  # Ensure it's a numpy array
    
    # Transform intensities to range from 0 to 1
    intensities_min = np.min(intensities)
    intensities_max = np.max(intensities)
    intensities = (intensities - intensities_min) / (intensities_max - intensities_min)
    
    # Mirror intensities across the y-axis
    intensities = 1 - intensities  # Flip the intensities
    
    # Plotting the spectrum
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wavenumbers, intensities, color='blue')
    ax.set_xlabel('Wavenumber (cm⁻¹)')
    ax.set_ylabel('Transmittance')
    ax.set_title(title)
    ax.invert_xaxis()
    ax.set_ylim(-0.05, 1.05)  # Add a 5% margin on top and bottom
    ax.grid(True)
    plt.show()

def identify_functional_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES string"

    functional_groups = []
    patterns = {
        "Alkyl halide": "[C][F,Cl,Br,I]",
        "Alcohol": "[OX2H]",
        "Aldehyde": "[CH1](=O)[#6]",
        "Ketone": "[#6][CX3](=O)[#6]",
        "Carboxylic acid": "[CX3](=O)[OX2H1]",
        "Ester": "[#6][CX3](=O)[OX2H0][#6]",
        "Ether": "[OD2]([#6])[#6]",
        "Amine": "[NX3;H2,H1,H0][#6]",
        "Amide": "[NX3][CX3](=[OX1])[#6]",
        "Nitro": "[N+](=O)[O-]",
        "Nitrile": "[C]#N",
        "Sulfide": "[#16X2H0]",
        "Sulfoxide": "[#16X3](=[OX1])",
        "Sulfone": "[#16X4](=[OX1])(=[OX1])",
        "Phosphate": "[P](=O)([O-])([O-])",
        "Phenol": "[OX2H][c]",
        "Imine": "[CX3]=[NX2]",
        "Alkene": "[CX2]=[CX2]",
        "Alkyne": "[CX2]#[CX2]",
        "Thiol": "[SX2H]",
        "Acyl chloride": "[CX3](=[OX1])[Cl]",
        "Anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
        "Lactam": "[NX3R][CX3R](=[OX1])[#6R]",
        "Aromatic": "[a]"
    }

    for name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            functional_groups.append(name)

    # Additional checks for specific structures
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[R]')) and not any(group in functional_groups for group in ["Aromatic", "Lactam"]):
        functional_groups.append("Alicyclic")

    return functional_groups

def generate_functional_group_vector(functional_groups):
    all_groups = [
        "Alkyl halide", "Alcohol", "Aldehyde", "Ketone", "Carboxylic acid",
        "Ester", "Ether", "Amine", "Amide", "Nitro", "Nitrile", "Sulfide", 
        "Sulfoxide", "Sulfone", "Phosphate", "Phenol", "Imine", "Alkene", 
        "Alkyne", "Thiol", "Acyl chloride", "Anhydride", "Lactam", 
        "Aromatic"
    ]
    
    vector = [1 if group in functional_groups else 0 for group in all_groups]
    return vector

def get_smiles_from_name(compound_name):
    # Try the original name first
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    # If not found, try removing parentheses and their contents
    modified_name = compound_name.split('(')[0].strip()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{modified_name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    # If still not found, try replacing spaces with hyphens
    modified_name = compound_name.replace(' ', '-')
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{modified_name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    return None

# Connect to the SQLite database
conn = sqlite3.connect('spectra.db')
cursor = conn.cursor()

# Fetch the first 100 molecules and their spectrum data
cursor.execute('''
    SELECT m.id, m.title, m.names, m.state, s.x, s.y
    FROM Molecules m
    JOIN SpectrumData s ON m.id = s.molecule_id
    WHERE m.id <= 100
    GROUP BY m.id
    ORDER BY m.id
''')

data = cursor.fetchall()

if data:
    all_datasets = []
    for molecule_id, title, names, state, x_values, y_values in data:
        # Fetch all x and y values for this molecule
        cursor.execute('''
            SELECT x, y
            FROM SpectrumData
            WHERE molecule_id = ?
            ORDER BY x
        ''', (molecule_id,))
        spectrum_data = cursor.fetchall()
        x_values, y_values = zip(*spectrum_data)

        # Create PyTorch tensors
        wavelengths_tensor = torch.tensor(x_values, dtype=torch.float32)
        intensities_tensor = torch.tensor(y_values, dtype=torch.float32)

        # Get SMILES from name and identify functional groups
        smiles = get_smiles_from_name(title)
        if smiles:
            functional_groups = identify_functional_groups(smiles)
            functional_group_vector = generate_functional_group_vector(functional_groups)
            functional_group_tensor = torch.tensor(functional_group_vector, dtype=torch.float32)
        else:
            print(f"Could not find SMILES for {names}")
            functional_group_tensor = torch.zeros(24, dtype=torch.float32)  # Assuming 24 functional groups

        class SpectrumDataset(torch.utils.data.Dataset):
            def __init__(self, wavelengths_tensor, intensities_tensor, functional_group_tensor):
                self.wavelengths_tensor = wavelengths_tensor
                self.intensities_tensor = intensities_tensor
                self.functional_group_tensor = functional_group_tensor

            def __len__(self):
                return len(self.wavelengths_tensor)

            def __getitem__(self, idx):
                return self.wavelengths_tensor[idx], self.intensities_tensor[idx], self.functional_group_tensor

        # Create the dataset
        dataset = SpectrumDataset(wavelengths_tensor, intensities_tensor, functional_group_tensor)
        all_datasets.append(dataset)

        print(f"Dataset created for molecule {molecule_id}: {names}")
        print(f"SMILES: {smiles}")
        print("Functional group vector:")
        print(dataset.functional_group_tensor)
        print("---")

    print(f"Total datasets created: {len(all_datasets)}")

    # Plot all spectra
    plt.figure(figsize=(15, 10))
    for dataset in all_datasets:
        plt.plot(dataset.wavelengths_tensor, dataset.intensities_tensor, alpha=0.5)
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('Intensity')
    plt.title('IR Spectra of First 100 Molecules')
    plt.show()

else:
    print("No data found in the database.")

# Close the database connection
conn.close()
