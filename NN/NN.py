import torch
import numpy as np
import matplotlib.pyplot as plt
import torch.nn.functional as F

data = torch.load('ir_spectroscopy_dataset.pt')
X_values = data['X_values']
Y_values = data['Y_values']
Labels = data['Labels']

print(X_values)


def preprocess_data(X_values, Y_values, Labels, cut_front=10, cut_end=10):
    """
    Preprocess the X_values and Y_values tensors by trimming a few points from the beginning and end,
    and removing data points with all NaN values or only NaN intensities.
    
    Parameters:
    - X_values: torch.Tensor
        The tensor containing the X values (wavenumbers) of the IR spectra
    - Y_values: torch.Tensor
        The tensor containing the Y values of the IR spectra
    - Labels: torch.Tensor
        The tensor containing the labels corresponding to Y_values
    - cut_front: int, default=10
        The number of data points to cut from the front of each spectrum
    - cut_end: int, default=10
        The number of data points to cut from the end of each spectrum
    
    Returns:
    - processed_X: torch.Tensor
        The processed X values tensor
    - processed_Y: torch.Tensor
        The processed Y values tensor
    - processed_Labels: torch.Tensor
        The processed labels tensor
    """
    # Trim the front and end of each spectrum
    processed_X = X_values[cut_front:-cut_end]
    processed_Y = Y_values[:, cut_front:-cut_end]
    
    # Create a mask for non-NaN rows and rows that are not all NaN
    mask = ~torch.isnan(processed_Y).all(dim=1) & ~torch.isnan(processed_Y).any(dim=1)
    
    # Apply the mask to Y_values and Labels
    processed_Y = processed_Y[mask]
    processed_Labels = Labels[mask]
    
    return processed_X, processed_Y, processed_Labels

# Example usage:
processed_X_values, processed_Y_values, processed_Labels = preprocess_data(X_values, Y_values, Labels)
print(f"Original X shape: {X_values.shape}, Processed X shape: {processed_X_values.shape}")
print(f"Original Y shape: {Y_values.shape}, Processed Y shape: {processed_Y_values.shape}")
print(f"Original labels: {Labels.shape}, Processed labels: {processed_Labels.shape}")

#----// NN //---------------------------


class NN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.g = torch.Generator().manual_seed(2147483647)

        self.W1 = torch.nn.Parameter(torch.randn((980, 2000), generator=self.g))
        self.b1 = torch.nn.Parameter(torch.randn(2000, generator=self.g))
        self.W2 = torch.nn.Parameter(torch.randn((2000, 500), generator=self.g))
        self.b2 = torch.nn.Parameter(torch.randn(500, generator=self.g))
        self.W3 = torch.nn.Parameter(torch.randn((500, 28), generator=self.g))
        self.b3 = torch.nn.Parameter(torch.randn(28, generator=self.g))

        print(f"The NN has: {sum(p.numel() for p in self.parameters())} parameters")

    def forward(self, Y_values, Labels):
        h = torch.tanh(Y_values @ self.W1 + self.b1)
        h2 = torch.tanh(h @ self.W2 + self.b2)
        logits = h2 @ self.W3 + self.b3
        out = torch.sigmoid(logits)  # Use sigmoid to get values between 0 and 1

        # Create a target tensor from Labels
        target = Labels.to(logits.device)

        # Compute binary cross-entropy loss
        loss = torch.nn.functional.binary_cross_entropy(out, target)

        return loss, out

    def __call__(self, data_tensor):
        h = torch.tanh(data_tensor @ self.W1 + self.b1)
        h2 = torch.tanh(h @ self.W2 + self.b2)
        logits = h2 @ self.W3 + self.b3
        out = torch.sigmoid(logits)  # Use sigmoid to get values between 0 and 1

        return out

nn_model = NN()
optimizer = torch.optim.SGD(nn_model.parameters(), lr=0.1)


# Check if CUDA is available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Move the model to the GPU
nn_model = nn_model.to(device)

# Move optimizer to GPU
optimizer = torch.optim.SGD(nn_model.parameters(), lr=0.1)

# Move data tensors to GPU
data_tensor = torch.tensor(processed_Y_values, dtype=torch.float32).to(device)
label_indices = torch.tensor(processed_Labels, dtype=torch.float32).to(device)

print("Model and data moved to GPU successfully.")



#----// NN //---------------------------



# Set up the training loop
data_tensor = torch.tensor(processed_Y_values, dtype=torch.float32)
label_indices = torch.tensor(processed_Labels, dtype=torch.float32)

epoch = 0
while True:
    optimizer.zero_grad()
    loss, _ = nn_model.forward(data_tensor, label_indices)
    loss.backward()
    
    print(f"Epoch {epoch+1}, Loss: {loss.item():.4f}")
    optimizer.step()
    
    if loss.item() < 1e-2:
        print("Reached Stop Point!!")
        print(f"Converged after {epoch+1} epochs.")
        break
    
    epoch += 1






#-------------------------------
"""
    Plot a single IR spectrum.
    
    Parameters:
    - X: array-like
        The wavenumbers (x-axis values)
    - Y: array-like
        The absorbance or transmittance values (y-axis values)
    - title: str, default='IR Spectrum'
        Title of the plot
    - xlabel: str, default='Wavenumber (cm⁻¹)'
        Label for the x-axis
    - ylabel: str, default='Absorbance'
        Label for the y-axis
    - filename: str, default='ir_spectrum.png'
        Filename to save the plot
"""
"""
def plot_ir(X, Y, title='IR Spectrum', xlabel='Wavenumber (cm⁻¹)', ylabel='Absorbance', filename='ir_spectrum.png'):
    plt.figure(figsize=(10, 6))
    plt.plot(X, Y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    
    print(f"IR spectrum plot saved as '{filename}'")

# Example usage:
plot_ir(processed_X_values, processed_Y_values[2])
"""
