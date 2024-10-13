import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np

# Directory containing the IR files
ir_directory = 'IR'

# Get a list of all files in the IR directory
ir_files = [f for f in os.listdir(ir_directory) if os.path.isfile(os.path.join(ir_directory, f))]

# Print the list of files (optional, for verification)
print("Files in the IR directory:")
for file in ir_files:
    #print(file)
    pass

print(ir_files[0])
# You can now use 'ir_files' list to process each file

with open('IR/' + ir_files[4], 'r') as file:
    # Read the first line
    first_line = file.readline().strip()
    print("First line:", first_line)
    # Reset file pointer to the beginning
    file.seek(0)
    
    # Search for NAME= or NAMES= or CAS NAME=
    for line in file:
        if line.startswith('##NAME=') or line.startswith('##NAMES=') or line.startswith('##CAS NAME='):
            # Extract and print the content after NAME=, NAMES=, or CAS NAME=
            name_content = line.split('=', 1)[1].strip()
            print(f"Found name: {name_content}")
            break
    else:
        print("No NAME=, NAMES=, or CAS NAME= field found in the file.")

    # Search for CAS REGISTRY NO=
    file.seek(0)
    for line in file:
        if line.startswith('##CAS REGISTRY NO='):
            cas_content = line.split('=', 1)[1].strip()
            print(f"Found CAS Registry Number: {cas_content}")
            break
    else:
        print("No CAS REGISTRY NO= field found in the file.")
    
    # Reset file pointer to the beginning
    file.seek(0)

    # Initialize lists to store x and y points
    x_points = []
    y_points = []

    # Flag to indicate when we've reached the data section
    data_section = False

    for line in file:
        if line.startswith('##XYDATA=(X++(Y..Y))'):
            data_section = True
            continue
        if data_section:
            # Check if the line starts with '##END='
            if line.startswith('##END='):
                break
            # Replace '-' with ' ' in the line
            line = line.replace('-', ' ')
            values = line.split()
            if values:
                x = float(values[0])
                for y in values[1:]:
                    x_points.append(x)
                    y_points.append(float(y))
                    x += 1  # Increment x for each y value

    print(f"Number of data points: {len(x_points)}")
    print("First few x points:", x_points[:5])
    print("First few y points:", y_points[:5])


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

    
    # Print some diagnostic information
    print(f"Wavelengths range: {np.min(wavelengths)} to {np.max(wavelengths)}")
    print(f"Intensities range: {np.min(intensities)} to {np.max(intensities)}")

plot_spectrum(x_points, y_points, title="Infrared Spectrum")


