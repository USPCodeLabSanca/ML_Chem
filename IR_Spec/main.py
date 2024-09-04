import requests
from nistchempy import Compound
import pandas as pd
import matplotlib.pyplot as plt

def get_ir_data(compound_name):
    # Search for the compound
    results = Compound.search(compound_name)
    
    if not results:
        print(f"No results found for {compound_name}")
        return None, None
    
    # Get the first result
    compound = results[0]
    
    # Get IR spectral data
    ir_data = compound.ir
    
    if not ir_data:
        print(f"No IR data available for {compound_name}")
        return None, None
    
    # Extract xy data
    x = ir_data['x']
    y = ir_data['y']
    
    # Create labels
    labels = {
        'name': compound.name,
        'formula': compound.formula,
        'cas': compound.cas_rn,
        'molecular_weight': compound.molecular_weight
    }
    
    return pd.DataFrame({'wavenumber': x, 'intensity': y}), labels

def plot_ir_spectrum(df, compound_name):
    plt.figure(figsize=(12, 6))
    plt.plot(df['wavenumber'], df['intensity'])
    plt.title(f"IR Spectrum of {compound_name}")
    plt.xlabel("Wavenumber (cm^-1)")
    plt.ylabel("Intensity")
    plt.gca().invert_xaxis()
    plt.show()

# Example usage
compound_name = "ethanol"
ir_data, labels = get_ir_data(compound_name)

if ir_data is not None:
    print("Labels:")
    for key, value in labels.items():
        print(f"{key}: {value}")
    
    print("\nIR Data (first 5 rows):")
    print(ir_data.head())
    
    # Plot the IR spectrum
    plot_ir_spectrum(ir_data, compound_name)
