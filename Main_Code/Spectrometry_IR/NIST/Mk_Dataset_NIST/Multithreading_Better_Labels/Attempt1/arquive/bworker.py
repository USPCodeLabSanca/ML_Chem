import sys
import time

import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem
import requests

import socks
import socket
import ssl



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
    conn = sqlite3.connect('final.db')
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
    db_path = os.path.join("Out", "labeled.db")
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
    #export_datasets(Dataset)

