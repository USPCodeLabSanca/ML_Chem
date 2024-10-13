import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from sklearn.preprocessing import StandardScaler
import torch

def parse_jdx(file_path):
    """Parse a JDX file and extract wavelengths and intensities."""
    wavelengths = []
    intensities = []
    with open(file_path, 'r') as file:
        data_started = False
        for line in file:
            if line.startswith('##XYDATA'):
                data_started = True
                continue
            if data_started:
                if line.strip().startswith('#'):
                    break
                parts = line.strip().split()
                wavelength = float(parts[0])
                for intensity in parts[1:]:
                    wavelengths.append(wavelength)
                    intensities.append(float(intensity))
                    wavelength += 0.935748  # Add DELTAX value
    return np.array(wavelengths), np.array(intensities)

def plot_spectrum(wavelengths, intensities):
    plt.figure(figsize=(12, 6))
    plt.plot(wavelengths, intensities)
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('Transmittance')
    plt.title('METHANE Infrared Spectrum')
    plt.tight_layout()
    plt.savefig('methane_spectrum.png')
    plt.close()

def save_spectrum(wavelengths, intensities, file_path):
    """Save the spectrum data as a PyTorch tensor in a .pt file."""
    spectrum_data = torch.tensor(np.column_stack((wavelengths, intensities)))
    torch.save(spectrum_data, file_path)
    print(f"Spectrum data saved as '{file_path}'")

def resample_and_trim(wavelengths, intensities, num_points=1000, trim=10):
    """Resample the spectrum to a specific number of points and trim the ends."""
    # Resample to 1000 points
    f = interpolate.interp1d(wavelengths, intensities)
    new_wavelengths = np.linspace(wavelengths.min(), wavelengths.max(), num_points)
    new_intensities = f(new_wavelengths)
    
    # Trim 10 points from beginning and end
    return new_wavelengths[trim:-trim], new_intensities[trim:-trim]

if __name__ == "__main__":
    jdx_file = 'example.jdx'
    wavelengths, intensities = parse_jdx(jdx_file)
    
    # Resample to 1000 points and trim
    wavelengths, intensities = resample_and_trim(wavelengths, intensities)
    
    plot_spectrum(wavelengths, intensities)
    print("Spectrum plot saved as 'methane_spectrum.png'")
    
    # Save the spectrum data as a .pt file
    save_spectrum(wavelengths, intensities, 'example.pt')

# To import the x and y data from the .pt file:
# loaded_data = torch.load('methane_spectrum.pt')
# x_data = loaded_data[:, 0]  # First column (wavelengths)
# y_data = loaded_data[:, 1]  # Second column (intensities)
