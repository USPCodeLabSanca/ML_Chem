import sqlite3
import matplotlib.pyplot as plt
import pandas as pd

def visualize_spectral_database():
    # Connect to the database
    conn = sqlite3.connect('spectral_database.db')
    
    # Get all table names from the database
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    
    for table in tables:
        table_name = table[0]
        
        # Query to get all data from the current table
        query = f"SELECT * FROM {table_name}"
        
        # Read the data into a pandas DataFrame
        df = pd.read_sql_query(query, conn)
        
        # Check if we have data
        if df.empty:
            print(f"No data found in the table {table_name}.")
            continue
        
        # Create a new figure for each table
        plt.figure(figsize=(12, 6))
        
        # Plot all columns as separate lines
        for column in df.columns:
            if df[column].dtype in ['float64', 'int64']:
                plt.plot(df.index, df[column], label=column)
        
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.title(f'Data from {table_name}')
        plt.legend()
        plt.grid(True)
        plt.show()
    
    # Close the database connection
    conn.close()

# Call the function to visualize the data
visualize_spectral_database()
