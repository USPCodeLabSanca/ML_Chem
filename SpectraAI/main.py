import sqlite3
from tabulate import tabulate

def visualize_db(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Get list of tables
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()

    for table in tables:
        table_name = table[0]
        print(f"\nTable: {table_name}")
        
        # Get table structure
        cursor.execute(f"PRAGMA table_info({table_name})")
        columns = cursor.fetchall()
        print("\nTable Structure:")
        print(tabulate(columns, headers=["ID", "Name", "Type", "NotNull", "DefaultValue", "PK"]))

        # Get row count
        cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
        row_count = cursor.fetchone()[0]
        print(f"\nTotal number of rows: {row_count}")

        # Get table data
        cursor.execute(f"SELECT * FROM {table_name} LIMIT 10")
        rows = cursor.fetchall()
        if rows:
            print("\nTable Data (up to 10 rows):")
            headers = [column[1] for column in columns]
            print(tabulate(rows, headers=headers))
        else:
            print("\nNo data in this table.")

    conn.close()

if __name__ == "__main__":
    db_path = "spectral_database.db"
    visualize_db(db_path)