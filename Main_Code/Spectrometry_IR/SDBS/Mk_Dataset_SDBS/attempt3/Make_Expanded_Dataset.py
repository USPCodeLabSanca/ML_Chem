import torch
import numpy as np
import random

def load_dataset(pt_path='ir_spectroscopy_dataset.pt'):
    """Load the existing dataset from a .pt file with safety considerations."""
    data = torch.load(pt_path, weights_only=True)  # Set weights_only as needed
    X = data['X_values']
    Y = data['Y_values']
    Labels = data['Labels']
    return X, Y, Labels

def add_gaussian_noise(Y, mean=0.0, std=0.01):
    """Add Gaussian noise to the spectra."""
    noise = torch.randn(Y.size()) * std + mean
    return Y + noise

def scale_spectra(Y, scale_range=(0.95, 1.05)):
    """Randomly scale the intensity of the entire spectrum."""
    scale = torch.empty(1).uniform_(*scale_range).item()  # Single scalar
    return Y * scale

def shift_wavelength(Y, shift_amount=5):
    """
    Shift the spectrum slightly by interpolating.
    Note: This is a simplistic approach and may require more sophisticated methods.
    """
    Y_np = Y.numpy()
    shifted_Y = np.roll(Y_np, shift_amount, axis=0)  # Adjust axis if necessary
    return torch.FloatTensor(shifted_Y)

def augment_data(Y, augment_funcs, num_augments=2):
    """
    Apply augmentation functions to each spectrum to generate new examples.
    
    Args:
        Y (torch.Tensor): Original spectra.
        augment_funcs (list): List of augmentation functions.
        num_augments (int): Number of augmented samples to generate per original sample.
    
    Returns:
        torch.Tensor: Augmented spectra.
    """
    augmented_Y = []
    for i in range(Y.size(0)):
        for _ in range(num_augments):
            augmented = Y[i].clone()  # Ensure original Y is not modified
            for func in augment_funcs:
                augmented = func(augmented)
            augmented_Y.append(augmented)
    return torch.stack(augmented_Y, dim=0)

def save_dataset(X, Y, Labels, save_path='expanded_ir_spectroscopy_dataset.pt'):
    """Save the augmented dataset to a .pt file."""
    torch.save({
        'X_values': X,
        'Y_values': Y,
        'Labels': Labels
    }, save_path)
    print(f"Expanded dataset saved as '{save_path}'")

def main():
    # Load original dataset
    X, Y, Labels = load_dataset()
    print(f"Original dataset: Y shape {Y.shape}, Labels shape {Labels.shape}")
    
    # Define augmentation functions
    augment_funcs = [add_gaussian_noise, scale_spectra]  # Exclude shift_wavelength if not needed
    
    # Generate augmented spectra
    augmented_Y = augment_data(Y, augment_funcs, num_augments=2)
    print(f"Augmented Y shape: {augmented_Y.shape}")
    
    # Duplicate labels for augmented data
    augmented_Labels = Labels.repeat(2, 1)
    
    # Concatenate original and augmented data
    combined_Y = torch.cat([Y, augmented_Y], dim=0)
    combined_Labels = torch.cat([Labels, augmented_Labels], dim=0)
    
    print(f"Combined dataset: Y shape {combined_Y.shape}, Labels shape {combined_Labels.shape}")
    
    # Validate dimensions
    assert combined_Y.size(0) == combined_Labels.size(0), "Mismatch between data and labels sizes."
    
    # Save the expanded dataset
    save_dataset(X, combined_Y, combined_Labels)

if __name__ == "__main__":
    main()