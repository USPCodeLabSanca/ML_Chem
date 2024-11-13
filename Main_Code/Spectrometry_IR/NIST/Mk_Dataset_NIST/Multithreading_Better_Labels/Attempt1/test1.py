import sqlite3
import sys


from pyCheckmol import CheckMol
cm = CheckMol()
import time

labels = []


def get_label(db_path, start_index, end_index):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Fetch SMILES formulas and molecule_id from the Molecules table within the given range
        cursor.execute('SELECT smiles, molecule_id FROM Molecules WHERE id BETWEEN ? AND ?', (start_index, end_index))
        molecules_data = cursor.fetchall()

        # Print each SMILES formula
        for molecule in molecules_data:
            print(molecule[0])
            try:
                smi = molecule[0]
                res = cm.functionalGroupSmiles(smiles=smi, isString=True, generate3D=False, justFGcode=False, returnDataframe=False,deleteTMP=False)
                print(res['Functional Group'])
                labels.append((molecule[1], res['Functional Group']))
                print((molecule[1], res['Functional Group']))
            except:
                pass

    except sqlite3.OperationalError as e:
        print(f"An error occurred: {e}")

    # Close the database connection
    conn.close()

# (6528, ['halogen deriv. '])

def write_labels(db_path, labels):
    # Connect to the SQLite database
    conn = sqlite3.connect(f'labeled_{slice_index}.db')
    cursor = conn.cursor()

    try:
        # Create the Labels table if it doesn't exist
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS Labels (
                molecule_id INTEGER PRIMARY KEY,
                label TEXT
            )
        ''')

        # Insert labels into the Labels table
        for molecule_id, label in labels:
            # Convert the list of labels to a comma-separated string
            label_str = ', '.join(label)
            cursor.execute('''
                INSERT OR REPLACE INTO Labels (molecule_id, label)
                VALUES (?, ?)
            ''', (molecule_id, label_str))

        # Commit the transaction
        conn.commit()
    except sqlite3.OperationalError as e:
        print(f"An error occurred: {e}")

    # Close the database connection
    conn.close()







if __name__ == '__main__':
    db_path = 'final.db'
    array = [int(arg) for arg in sys.argv[1:]]
    start_index = array[0]
    end_index = array[1]
    #end_index = array[0] + 10

    slice_index = array[2]

    print(f"Process started with array: {array}")
    get_label(db_path, start_index, end_index)

    print("--------------// Final Steps //---,----------------")
    write_labels(db_path, labels)
