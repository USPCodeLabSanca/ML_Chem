import os
import re
import sqlite3

# Directory containing the .jdx files
dir_path = 'IR'  # Adjusted directory path

# Regular expression patterns
number_pattern_str = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
number_pattern = re.compile(r'^' + number_pattern_str + r'$')
range_pattern = re.compile(r'^(' + number_pattern_str + r')-(' + number_pattern_str + r')$')

# Connect to SQLite database (or create it if it doesn't exist)
conn = sqlite3.connect('spectra.db')
cursor = conn.cursor()

# Create tables
cursor.execute('''
    CREATE TABLE IF NOT EXISTS Molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        title TEXT,
        names TEXT,
        state TEXT
    )
''')

cursor.execute('''
    CREATE TABLE IF NOT EXISTS SpectrumData (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INTEGER,
        x REAL,
        y REAL,
        FOREIGN KEY (molecule_id) REFERENCES Molecules(id)
    )
''')

# Traverse the IR directory
for filename in os.listdir(dir_path):
    if filename.endswith('.jdx') or filename.endswith('.jsx'):
        filepath = os.path.join(dir_path, filename)
        print(f"Processing file: {filename}")
        title = ''
        names = ''
        state = ''
        xydata = []
        reading_xydata = False
        reading_peak_table = False
        delta_x = None

        # Open the file with the correct encoding
        with open(filepath, 'r', encoding='latin-1') as f:
            lines = f.readlines()

        for line in lines:
            line = line.strip()
            if line.startswith('##TITLE='):
                title = line[len('##TITLE='):].strip()
            elif line.startswith('##NAMES='):
                names = line[len('##NAMES='):].strip()
            elif line.startswith('##STATE='):
                state = line[len('##STATE='):].strip()
            elif line.startswith('##DELTAX='):
                delta_x_str = line[len('##DELTAX='):].strip()
                try:
                    delta_x = float(delta_x_str)
                except ValueError:
                    print(f"Warning: DELTAX value '{delta_x_str}' is not a valid float in file {filename}")
                    delta_x = None
            elif line.startswith('##XYDATA='):
                xydata_header = line[len('##XYDATA='):].strip()
                reading_xydata = True
                reading_peak_table = False
                xydata = []
            elif line.startswith('##PEAK TABLE='):
                peak_table_header = line[len('##PEAK TABLE='):].strip()
                reading_peak_table = True
                reading_xydata = False
                xydata = []
            elif line.startswith('##END='):
                reading_xydata = False
                reading_peak_table = False
            elif reading_xydata or reading_peak_table:
                if line == '':
                    continue
                tokens = line.strip().split()
                if reading_xydata:
                    if len(tokens) >= 2:
                        x_value_str = tokens[0]
                        # Handle x-value ranges
                        match = range_pattern.match(x_value_str)
                        if match:
                            x_start_str, x_end_str = match.group(1), match.group(2)
                            try:
                                x_start = float(x_start_str)
                                x_end = float(x_end_str)
                                x_value = (x_start + x_end) / 2
                            except ValueError:
                                print(f"Warning: Cannot convert x range values to float in file {filename}: '{x_start_str}', '{x_end_str}'")
                                continue
                        elif number_pattern.match(x_value_str):
                            x_value = float(x_value_str)
                        else:
                            print(f"Warning: Cannot parse x_value '{x_value_str}' in file {filename}")
                            continue
                        y_values = []
                        for y_str in tokens[1:]:
                            y_str_clean = y_str.strip()
                            # Handle y-value ranges
                            match = range_pattern.match(y_str_clean)
                            if match:
                                y_start_str, y_end_str = match.group(1), match.group(2)
                                try:
                                    y_start = float(y_start_str)
                                    y_end = float(y_end_str)
                                    y_value = (y_start + y_end) / 2
                                except ValueError:
                                    print(f"Warning: Cannot convert y range values to float in file {filename}: '{y_start_str}', '{y_end_str}'")
                                    continue
                            elif number_pattern.match(y_str_clean):
                                y_value = float(y_str_clean)
                            else:
                                print(f"Warning: Cannot parse y_value '{y_str}' in file {filename}")
                                continue
                            y_values.append(y_value)
                        x = x_value
                        for y in y_values:
                            xydata.append((x, y))
                            if delta_x is not None:
                                x += delta_x
                    else:
                        print(f"Warning: unexpected format in data line in file {filename}: {line}")
                elif reading_peak_table:
                    # In peak tables, data might be space-separated x,y pairs
                    for token in tokens:
                        # Handle x,y pairs
                        if ',' in token:
                            x_str, y_str = token.split(',', 1)
                            x_str_clean = x_str.strip()
                            y_str_clean = y_str.strip()
                            # Handle x-value ranges
                            match_x = range_pattern.match(x_str_clean)
                            if match_x:
                                x_start_str, x_end_str = match_x.group(1), match_x.group(2)
                                try:
                                    x_start = float(x_start_str)
                                    x_end = float(x_end_str)
                                    x_value = (x_start + x_end) / 2
                                except ValueError:
                                    print(f"Warning: Cannot convert x range values to float in file {filename}: '{x_start_str}', '{x_end_str}'")
                                    continue
                            elif number_pattern.match(x_str_clean):
                                x_value = float(x_str_clean)
                            else:
                                print(f"Warning: Cannot parse x_value '{x_str}' in file {filename}")
                                continue
                            # Handle y-value ranges
                            match_y = range_pattern.match(y_str_clean)
                            if match_y:
                                y_start_str, y_end_str = match_y.group(1), match_y.group(2)
                                try:
                                    y_start = float(y_start_str)
                                    y_end = float(y_end_str)
                                    y_value = (y_start + y_end) / 2
                                except ValueError:
                                    print(f"Warning: Cannot convert y range values to float in file {filename}: '{y_start_str}', '{y_end_str}'")
                                    continue
                            elif number_pattern.match(y_str_clean):
                                y_value = float(y_str_clean)
                            else:
                                print(f"Warning: Cannot parse y_value '{y_str}' in file {filename}")
                                continue
                            xydata.append((x_value, y_value))
                        else:
                            # Single x-value, y defaults to 1.0
                            x_str_clean = token.strip()
                            if number_pattern.match(x_str_clean):
                                x_value = float(x_str_clean)
                                y_value = 1.0
                                xydata.append((x_value, y_value))
                            else:
                                print(f"Warning: Cannot parse peak value '{token}' in file {filename}")
            # Else: other lines, ignore

        if state.lower() != 'gas':
            print(f"Skipping molecule '{title}' as it's not in gas state.")
            continue
        if title == '' and names == '':
            print(f"Warning: No title or names found in file {filename}")
            continue
        print(f"Processed data for molecule '{title}'.")
        print(f"State: {state}")

        # Insert molecule data into the database
        cursor.execute('''
            INSERT INTO Molecules (title, names, state)
            VALUES (?, ?, ?)
        ''', (title, names, state))
        molecule_id = cursor.lastrowid

        # Insert spectrum data into the database
        spectrum_data = [(molecule_id, x, y) for x, y in xydata]
        cursor.executemany('''
            INSERT INTO SpectrumData (molecule_id, x, y)
            VALUES (?, ?, ?)
        ''', spectrum_data)

        # Commit after each molecule to save data incrementally
        conn.commit()

# Close the database connection
conn.close()
print("Data has been successfully stored in 'spectra.db'.")
