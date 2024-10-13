import sqlite3
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import time
import json

def create_labeled_database():
    with sqlite3.connect('labeled.db', timeout=10) as conn:
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS MoleculeData
                     (id INTEGER PRIMARY KEY, name TEXT, vector TEXT)''')
        c.execute('''CREATE TABLE IF NOT EXISTS SpectrumData
                     (id INTEGER PRIMARY KEY, molecule_id INTEGER, 
                      wavelength REAL, intensity REAL,
                      FOREIGN KEY(molecule_id) REFERENCES MoleculeData(id))''')

def get_smiles_from_name(name):
    # Try the original name first
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    # If not found, try removing parentheses and their contents
    modified_name = name.split('(')[0].strip()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{modified_name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    # If still not found, try replacing spaces with hyphens
    modified_name = name.replace(' ', '-')
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{modified_name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    return None

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

def analyze_molecule_and_store(name):
    smiles = get_smiles_from_name(name)
    if smiles:
        functional_groups = identify_functional_groups(smiles)
        vector = generate_functional_group_vector(functional_groups)
        vector_str = ','.join(map(str, vector))
        
        max_retries = 5
        for attempt in range(max_retries):
            try:
                with sqlite3.connect('labeled.db', timeout=10) as conn:
                    c = conn.cursor()
                    c.execute("INSERT INTO MoleculeData (name, vector) VALUES (?, ?)", (name, vector_str))
                break
            except sqlite3.OperationalError as e:
                if attempt < max_retries - 1:
                    time.sleep(1)
                else:
                    print(f"Failed to insert data for {name} after {max_retries} attempts: {e}")
                    return False
        
        print(f"Molecule: {name}")
        print(f"SMILES: {smiles}")
        print("Functional groups:")
        for group in functional_groups:
            print(f"- {group}")
        print(f"Vector: {vector}")
        print()
        return True
    else:
        print(f"Could not find SMILES for {name}")
        print()
        return False

def transfer_data_from_spectral_to_labeled():
    with sqlite3.connect('spectral_database.db', timeout=10) as spectral_conn:
        spectral_cursor = spectral_conn.cursor()
        spectral_cursor.execute("SELECT id, name FROM SpectraNames")
        molecules = spectral_cursor.fetchall()
    
    for molecule_id, molecule_name in molecules:
        if not analyze_molecule_and_store(molecule_name):
            print(f"Skipping spectrum data transfer for {molecule_name} due to analysis failure.")
            continue
        
        with sqlite3.connect('spectral_database.db', timeout=10) as spectral_conn:
            spectral_cursor = spectral_conn.cursor()
            spectral_cursor.execute("SELECT wavelength, intensity FROM SpectraData WHERE name_id = ?", (molecule_id,))
            spectra_data = spectral_cursor.fetchall()
        
        max_retries = 5
        for attempt in range(max_retries):
            try:
                with sqlite3.connect('labeled.db', timeout=10) as labeled_conn:
                    labeled_cursor = labeled_conn.cursor()
                    labeled_cursor.execute("SELECT id FROM MoleculeData WHERE name = ?", (molecule_name,))
                    result = labeled_cursor.fetchone()
                    
                    if result is None:
                        print(f"No data found in labeled database for {molecule_name}. Skipping spectrum data transfer.")
                        break
                    
                    molecule_data_id = result[0]
                    
                    labeled_cursor.executemany("INSERT INTO SpectrumData (molecule_id, wavelength, intensity) VALUES (?, ?, ?)",
                                               [(molecule_data_id, wavelength, intensity) for wavelength, intensity in spectra_data])
                break
            except sqlite3.OperationalError as e:
                if attempt < max_retries - 1:
                    time.sleep(1)
                else:
                    print(f"Failed to update data for {molecule_name} after {max_retries} attempts: {e}")

if __name__ == "__main__":
    create_labeled_database()
    transfer_data_from_spectral_to_labeled()
