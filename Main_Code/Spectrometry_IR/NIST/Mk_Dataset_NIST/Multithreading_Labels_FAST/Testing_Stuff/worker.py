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


def process_array(array):
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

    else:
        print("No data found in the database for the given range.")

    # Close the database connection
    conn.close()

    return all_datasets


def export_datasets(all_datasets, slice_index):
    combined_dataset = torch.utils.data.ConcatDataset(all_datasets)
    
    # Create a dictionary to store all the data
    data_dict = {
        'wavelengths': [],
        'intensities': [],
        'functional_groups': []
    }
    
    # Iterate through the combined dataset and collect all data
    for i in range(len(combined_dataset)):
        wavelengths, intensities, functional_groups = combined_dataset[i]
        data_dict['wavelengths'].append(wavelengths)
        data_dict['intensities'].append(intensities)
        data_dict['functional_groups'].append(functional_groups)
    
    # Convert lists to tensors
    for key in data_dict:
        data_dict[key] = torch.stack(data_dict[key])
    
    # Save the dictionary as a .pt file
    filename = f"dataset_slice__{slice_index}.pt"
    torch.save(data_dict, filename)
    print(f"Dataset exported to {filename}")







if __name__ == '__main__':
    # Skip the script name (sys.argv[0])
    array = [int(arg) for arg in sys.argv[1:]]
    
    print(f"Process started with array: {array}")
    # Process the array
    Dataset = process_array(array)
    export_datasets(Dataset, array[2])

    print(f"Process completed. Input: {array}, Output: {len(Dataset.wavelengths)}")
