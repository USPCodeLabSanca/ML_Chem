import sqlite3
import glob
import sys

def combine_label_databases():
    try:
        # Create a new database for the combined labels
        conn_final = sqlite3.connect('labeled.db')
        cursor_final = conn_final.cursor()

        # Create the Labels table in the final database
        cursor_final.execute('''
            CREATE TABLE IF NOT EXISTS Labels (
                molecule_id INTEGER PRIMARY KEY,
                label TEXT
            )
        ''')

        # Get list of all labeled_*.db files, excluding 'labeled.db'
        labeled_dbs = [db for db in glob.glob('labeled_*.db') if db != 'labeled.db']

        print(f"Found {len(labeled_dbs)} labeled database files.")
        if len(labeled_dbs) == 0:
            print("No source databases found. Exiting.")
            return

        print("Source databases to be processed:")
        for db in labeled_dbs:
            print(f"- {db}")

        # Dictionary to store all records
        all_records = {}

        # Combine data from all labeled_*.db files
        for db_file in labeled_dbs:
            print(f"\nProcessing {db_file}")
            conn = sqlite3.connect(db_file)
            cursor = conn.cursor()

            # Verify the number of records
            cursor.execute('SELECT COUNT(*) FROM Labels')
            count = cursor.fetchone()[0]
            print(f"Number of records in {db_file}: {count}")

            # Get all records from current db
            cursor.execute('SELECT molecule_id, label FROM Labels')
            records = cursor.fetchall()

            print(f"Retrieved {len(records)} records from {db_file}")

            # Store records in dictionary, preserving all labels for each molecule_id
            for molecule_id, label in records:
                if molecule_id in all_records:
                    # Combine labels if molecule_id already exists
                    existing_labels = set(all_records[molecule_id].split(', '))
                    new_labels = set(label.split(', '))
                    combined_labels = existing_labels.union(new_labels)
                    all_records[molecule_id] = ', '.join(sorted(combined_labels))
                else:
                    all_records[molecule_id] = label

            conn.close()

            # Debugging: After each file
            print(f"Current total unique molecule_ids after {db_file}: {len(all_records)}")

        print(f"\nTotal unique molecule_ids collected: {len(all_records)}")

        if len(all_records) == 0:
            print("No records to insert into labeled.db. Exiting.")
            return

        # Insert combined records into final db in batches
        batch_size = 1000
        records_items = list(all_records.items())
        total_records = len(records_items)

        print("Starting batch insertion into labeled.db...")
        for start in range(0, total_records, batch_size):
            batch = records_items[start:start + batch_size]
            cursor_final.executemany('''
                INSERT OR REPLACE INTO Labels (molecule_id, label)
                VALUES (?, ?)
            ''', batch)
            conn_final.commit()
            print(f"Inserted records {start + 1} to {min(start + batch_size, total_records)}")

        # Close final database
        conn_final.close()

        print("\nAll label databases have been successfully combined into labeled.db")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    combine_label_databases()
