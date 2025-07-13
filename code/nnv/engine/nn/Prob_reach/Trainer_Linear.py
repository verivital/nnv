import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import scipy.io
import numpy as np
import random
import time
import gc


# Define the path to the .mat file
mat_file_path = r"Temp_files_mid_run/Reduced_dimension.mat"

# Load MATLAB data
mat_data = scipy.io.loadmat(mat_file_path)

# Extract and transpose data for PyTorch
x = torch.tensor(mat_data['X'].T, dtype=torch.float32)  # Shape [10000, 5376]
y = torch.tensor(mat_data['dYV'].T, dtype=torch.float32)  # Shape [10000, 10]
epochs = int(mat_data['epochs'].flatten()[0])
lr = mat_data['lr'].flatten()[0]





# Compute mean and std for normalization
y_mean = y.mean(dim=0, keepdim=True)  # Shape [1, 10]
y_std = y.std(dim=0, keepdim=True) + 1e-8  # Avoid division by zero

# Normalize y
y_norm = (y - y_mean) / y_std  # Shape [10000, 10]

# Check for GPU availability
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Create DataLoader for mini-batch training
batch_size = 20
dataset = TensorDataset(x, y_norm)
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)


class LinearModel(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(LinearModel, self).__init__()
        self.linear = nn.Linear(input_dim, output_dim)  # y = W*x + b

    def forward(self, x):
        return self.linear(x)

# Define Model and move to GPU

input_dim   = x.shape[1]

output_dim  = y.shape[1]

model = LinearModel(input_dim, output_dim).to(device)

# Define Loss and Optimizer
criterion = nn.MSELoss(reduction='sum')
optimizer = optim.Adam(model.parameters(), lr=lr)


start_time =time.time()

for epoch in range(epochs):
    total_loss = 0
    # uuj = 0
    for x_batch, y_batch in dataloader:
        # uuj = uuj+1
        # print(uuj)
        x_batch, y_batch = x_batch.to(device), y_batch.to(device)

        optimizer.zero_grad()
        y_pred = model(x_batch)
        loss = criterion(y_pred, y_batch)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    if epoch % 10 == 0:
        print(f'Epoch [{epoch}/{epochs}], Loss: {total_loss / len(dataloader):.4f}')
        


# Extract and Save Optimized Weights and Biases
W = model.linear.weight.detach().cpu().numpy()  # Shape [10, 5376]
b = model.linear.bias.detach().cpu().numpy()  # Shape [10]

# Reshape y_std to [10, 1] for correct broadcasting
y_std_np = y_std.numpy().reshape(-1, 1)

# Denormalize final layer (W3, b3)
W_denorm = (y_std_np * W)
b_denorm = (y_std.numpy() * b) + y_mean.numpy()


end_time = time.time()

training_time = end_time - start_time

# Save to MATLAB format
save_path = r"Temp_files_mid_run/trained_Linear_weights_norm.mat"
scipy.io.savemat(save_path, {
    'W': W_denorm, 'b': b_denorm, 'Model_training_time': training_time
})

print(f"Trained parameters saved to {save_path}")


del x, y, x_batch, y_batch, y_pred 
gc.collect()
torch.cuda.empty_cache()