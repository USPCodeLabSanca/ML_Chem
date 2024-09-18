import requests
from rdkit import Chem
from rdkit.Chem import Descriptors

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

    # Define SMARTS patterns for various functional groups
    patterns = {
        "Alkyl halide": "[C][Cl,Br,I,F]",
        "Alcohol": "[OX2H]",
        "Aldehyde": "[CX3H1](=O)[#6]",
        "Ketone": "[#6][CX3](=O)[#6]",
        "Carboxylic acid": "[CX3](=O)[OX2H1]",
        "Ester": "[#6][CX3](=O)[OX2H0][#6]",
        "Ether": "[OD2]([#6])[#6]",
        "Amine": "[NX3;H2,H1;!$(NC=O)]",
        "Amide": "[NX3][CX3](=[OX1])[#6]",
        "Nitro": "[N+](=O)[O-]",
        "Nitrile": "[C]#N",
        "Sulfide": "[#16X2H0]",
        "Sulfoxide": "[#16X3](=[OX1])",
        "Sulfone": "[#16X4](=[OX1])(=[OX1])",
        "Phosphate": "[P](=O)([O-])([O-])",
        "Phenol": "[OX2H][cX3]:[c]",
        "Imine": "[CX3]=[NX2]",  # Added Imine pattern
        "Alkene": "[CX2]=[CX2]",  # Added Alkene pattern
        "Alkyne": "[CX2]#[CX2]",  # Added Alkyne pattern
        "Thiol": "[SX2H]",  # Added Thiol pattern
        "Acyl chloride": "[CX3](=[OX1])[Cl]",  # Added Acyl chloride pattern
        "Anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",  # Added Anhydride pattern
        "Epoxide": "[OX2r3]1[CX4r3][CX4r3]1",  # Added Epoxide pattern
    }

    for name, smarts in patterns.items():
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                functional_groups.append(name)
        except Exception as e:
            print(f"Error processing pattern for {name}: {e}")

    return functional_groups

def analyze_molecule(name):
    smiles = get_smiles_from_name(name)
    if smiles:
        functional_groups = identify_functional_groups(smiles)
        print(f"Molecule: {name}")
        print(f"SMILES: {smiles}")
        print("Functional groups:")
        for group in functional_groups:
            print(f"- {group}")
        print()
    else:
        print(f"Could not find SMILES for {name}")
        print()

if __name__ == "__main__":
    molecules = [
        "1,2-Dichloropropane",
        "1,3-Dichloropropene(cis)",
        "1-Naphthol",
        "1-Naphthyl-acetic-acid",
        "2,4'-DDT",
        "2,4,5-Trichlorophenol",
        "2,4-DB",
        "2,6-Dichlorobenzamide",
        "4,4'-DDT",
        "Acephate",
        "Alachlor",
        "Aldicarb",
        "Aldrin",
        "Alloxydim-sodium",
        "Amitraz"
    ]

    for molecule in molecules:
        analyze_molecule(molecule)