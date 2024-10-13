import sys
import time

import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem
import requests




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






class SpectrumDataset(torch.utils.data.Dataset):
    def __init__(self):
        self.ids = []
        self.smiles = []
        self.wavelengths_tensor = []
        self.intensities_tensor = []
        self.functional_group_tensor = []
        self.names = []

    def __len__(self):
        return len(self.wavelengths_tensor)

    def __getitem__(self, idx):
        return self.wavelengths_tensor[idx], self.intensities_tensor[idx], self.functional_group_tensor[idx], self.names[idx]

    def add_item(self, wavelengths, intensities, functional_group, name, smiles, id):
        self.wavelengths_tensor.append(wavelengths)
        self.intensities_tensor.append(intensities)
        self.functional_group_tensor.append(functional_group)
        self.names.append(name)
        self.smiles.append(smiles)
        self.ids.append(id)
        return self


def process_array(Dataset, array):
    start_index, end_index, slice_index = array
    # Connect to the SQLite database
    conn = sqlite3.connect('spectra.db')
    cursor = conn.cursor()

    # Fetch molecules and their spectrum data within the given range
    cursor.execute('''
        SELECT m.id, m.title, m.names, m.state, s.x, s.y
        FROM Molecules m
        JOIN SpectrumData s ON m.id = s.molecule_id
        WHERE m.id BETWEEN ? AND ?
        GROUP BY m.id
        ORDER BY m.id
    ''', (start_index, end_index))

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


            Dataset = Dataset.add_item(wavelengths_tensor, intensities_tensor, functional_group_tensor, title, smiles, molecule_id)

            print(f"Dataset created for molecule {molecule_id}: {names}")
            print(f"SMILES: {smiles}")
            print("Functional group vector:")
            print(Dataset.functional_group_tensor[-1])
            print("---")

        print(f"Total datasets created: {len(Dataset.wavelengths_tensor)}")

    else:
        print("No data found in the database for the given range.")

    # Close the database connection
    conn.close()

    return Dataset




def export_datasets(Dataset):
    import sqlite3
    import os

    # Connect to the SQLite database
    db_path = os.path.join("Out", "final.db")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create the tables if they don't exist
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS Molecules (
        id INTEGER PRIMARY KEY,
        molecule_id INTEGER,
        name TEXT,
        smiles TEXT
    )
    ''')

    cursor.execute('''
    CREATE TABLE IF NOT EXISTS SpectrumData (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INTEGER,
        wavelength REAL,
        intensity REAL,
        FOREIGN KEY (molecule_id) REFERENCES Molecules(id)
    )
    ''')

    cursor.execute('''
    CREATE TABLE IF NOT EXISTS FunctionalGroups (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INTEGER,
        group_name TEXT,
        present INTEGER,
        FOREIGN KEY (molecule_id) REFERENCES Molecules(id)
    )
    ''')

    # Prepare and insert the data
    for molecule_id in Dataset.ids:
        index = Dataset.ids.index(molecule_id)
        # Check if molecule already exists
        cursor.execute('SELECT id FROM Molecules WHERE molecule_id = ?', (molecule_id,))
        existing_molecule = cursor.fetchone()

        if not existing_molecule:
            # Insert molecule data
            cursor.execute('''
            INSERT OR IGNORE INTO Molecules (molecule_id, name, smiles)
            VALUES (?, ?, ?)
            ''', (molecule_id, Dataset.names[index], Dataset.smiles[index]))
            
            # Insert spectrum data
            spectrum_data = [(molecule_id, float(w), float(i)) for w, i in zip(Dataset.wavelengths_tensor[index], Dataset.intensities_tensor[index])]
            cursor.executemany('''
            INSERT OR IGNORE INTO SpectrumData (molecule_id, wavelength, intensity)
            VALUES (?, ?, ?)
            ''', spectrum_data)

            # Insert functional group data
            functional_groups = [
                "Alkyl halide", "Alcohol", "Aldehyde", "Ketone", "Carboxylic acid",
                "Ester", "Ether", "Amine", "Amide", "Nitro", "Nitrile", "Sulfide", 
                "Sulfoxide", "Sulfone", "Phosphate", "Phenol", "Imine", "Alkene", 
                "Alkyne", "Thiol", "Acyl chloride", "Anhydride", "Lactam", 
                "Aromatic"
            ]
            functional_group_data = [(molecule_id, group, int(present)) for group, present in zip(functional_groups, Dataset.functional_group_tensor[index])]
            cursor.executemany('''
            INSERT OR IGNORE INTO FunctionalGroups (molecule_id, group_name, present)
            VALUES (?, ?, ?)
            ''', functional_group_data)

    # Commit the changes and close the connection
    conn.commit()
    conn.close()

    print(f"Dataset exported to {db_path}")



if __name__ == '__main__':
    array = [int(arg) for arg in sys.argv[1:]]
    slice_index = array[2]

    Dataset = SpectrumDataset()
    print(f"Process started with array: {array}")
    # Process the array
    Dataset = process_array(Dataset, array)
    # Processing the data

    print("--------------// Final Steps //---,----------------")
    export_datasets(Dataset)

