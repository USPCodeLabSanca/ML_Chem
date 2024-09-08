import os
import pandas as pd
import matplotlib.pyplot as plt

IR_SPECTRA_DIR = "ir_spectra"

def plot_ir_spectrum(cas, x_values, y_values, output_dir=IR_SPECTRA_DIR):
    """Plot the IR spectrum for a given compound."""
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values)
    plt.title(f"IR Spectrum for CAS {cas}")
    plt.xlabel("Wavenumber (cm⁻¹)")
    plt.ylabel("Transmittance (%)")
    plt.savefig(f"{output_dir}/{cas}_ir_spectrum.png")
    plt.close()

def create_ir_dataset(compound_data_file, output_dir=IR_SPECTRA_DIR):
    """Create a dataset of IR spectra from local CSV files."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Read the compound data
    compounds_df = pd.read_csv(compound_data_file)
    
    for _, compound in compounds_df.iterrows():
        cas = compound['CAS']
        ir_file = f"{IR_SPECTRA_DIR}/{cas}_ir_spectrum.csv"
        
        if os.path.exists(ir_file):
            ir_data = pd.read_csv(ir_file)
            plot_ir_spectrum(cas, ir_data['X'], ir_data['Y'], output_dir)
        else:
            print(f"IR spectrum file not found for CAS {cas}")

    print(f"IR spectra plots saved in {output_dir} directory")

if __name__ == '__main__':
    compound_data_file = "compound_dataset.csv"
    create_ir_dataset(compound_data_file)
import requests
import pandas as pd

def get_ir_spectroscopy_data(compound_id):
    # API endpoint (replace with actual database URL)
    url = f"https://api.example.com/ir_spectroscopy/{compound_id}"
    
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        
        # Convert to DataFrame
        df = pd.DataFrame(data, columns=['Wavenumber', 'Transmittance'])
        return df
    else:
        print(f"Error: Unable to retrieve data (Status code: {response.status_code})")
        return None

# Example usage
compound_id = "C6H6"  # Benzene
ir_data = get_ir_spectroscopy_data(compound_id)

if ir_data is not None:
    print(ir_data.head())
    # Further processing or visualization can be done here