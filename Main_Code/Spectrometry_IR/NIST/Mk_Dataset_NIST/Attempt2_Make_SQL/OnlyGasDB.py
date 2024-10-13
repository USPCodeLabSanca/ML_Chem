import os
import re
import matplotlib.pyplot as plt
import numpy as np

# Directory containing the .jdx files
dir_path = 'IR'  # Adjusted directory path

# Regular expression patterns
number_pattern_str = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
number_pattern = re.compile(r'^' + number_pattern_str + r'$')
range_pattern = re.compile(r'^(' + number_pattern_str + r')-(' + number_pattern_str + r')$')

def plot_spectrum(wavelengths, intensities, title="Infrared Spectrum"):
    wavelengths = np.array(wavelengths)
    wavenumbers = wavelengths
    intensities = np.array(intensities)  # Ensure it's a numpy array
    
    # Transform intensities to range from 0 to 1
    intensities_min = np.min(intensities)
    intensities_max = np.max(intensities)
    intensities = (intensities - intensities_min) / (intensities_max - intensities_min)
    
    # Mirror intensities across the y-axis
    intensities = 1 - intensities  # Flip the intensities
    
    # Plotting the spectrum
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wavenumbers, intensities, color='blue')
    ax.set_xlabel('Wavenumber (cm⁻¹)')
    ax.set_ylabel('Transmittance')
    ax.set_title(title)
    ax.invert_xaxis()
    ax.set_ylim(-0.05, 1.05)  # Add a 5% margin on top and bottom
    ax.grid(True)
    plt.show()

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
        with open(filepath, 'r') as f:
            lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith('##TITLE='):
                title = line[len('##TITLE='):].strip()
            elif line.startswith('##NAMES='):
                names = line[len('##NAMES='):].strip()
            elif line.startswith('##STATE='):
                state = line[len('##STATE='):].strip()
                if state.lower() != 'gas':
                    break  # Skip this molecule if it's not in gas state
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

        # Plot the spectrum
        wavelengths = np.array([x for x, y in xydata])
        intensities = np.array([y for x, y in xydata])
        plot_spectrum(wavelengths, intensities, title=title)
