import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import scipy.io
import numpy as np
import random
import time
import gc

def estimate_lipschitz(x, y, num_samples=1000):
    n = x.shape[0]
    print(n)
    slopes = []

    for _ in range(num_samples):
        i, j = random.sample(range(n), 2)
        diff_x = x[i] - x[j]
        diff_y = y[i] - y[j]

        norm_x = torch.norm(diff_x)
        norm_y = torch.norm(diff_y)

        if norm_x > 1e-8:  # avoid division by near-zero
            slope = norm_y / norm_x
            slopes.append(slope.item())

    return max(slopes)


# Define the path to the .mat file
mat_file_path = r"Temp_files_mid_run/Reduced_dimension.mat"

# Load MATLAB data
mat_data = scipy.io.loadmat(mat_file_path)

# Extract and transpose data for PyTorch
x = torch.tensor(mat_data['X'].T, dtype=torch.float32)  # Shape [10000, 5376]
y = torch.tensor(mat_data['dYV'].T, dtype=torch.float32)  # Shape [10000, 10]
dims = mat_data['dims'].flatten().astype(int).tolist()
epochs = int(mat_data['epochs'].flatten()[0])

# Estimate λ before training
lam = max( 10.0 , 5*estimate_lipschitz(x, y) )
print(f"Estimated Lipschitz constant (empirical): {lam:.4f}")





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

# Define Neural Network Model with Two Hidden Layers
class ReLUNetwork(nn.Module):
    def __init__(self, input_dim, hidden_dim1, hidden_dim2, output_dim):
        super(ReLUNetwork, self).__init__()
        self.hidden1 = nn.Linear(input_dim, hidden_dim1)
        self.relu = nn.ReLU()
        self.hidden2 = nn.Linear(hidden_dim1, hidden_dim2)
        self.output = nn.Linear(hidden_dim2, output_dim)

    def forward(self, x):
        x = self.hidden1(x)
        x = self.relu(x)
        x = self.hidden2(x)
        x = self.relu(x)
        x = self.output(x)
        return x

# Define Model and move to GPU

input_dim   = x.shape[1]
if dims[0] == -1:
    hidden_dim1 = input_dim
else:
    hidden_dim1 = dims[0]

output_dim  = y.shape[1]
if dims[1] == -1:
    hidden_dim2 = output_dim
else:
    hidden_dim2 = dims[1]

model = ReLUNetwork(input_dim, hidden_dim1, hidden_dim2, output_dim).to(device)

# Define Loss and Optimizer
criterion = nn.MSELoss(reduction='sum')
optimizer = optim.Adam(model.parameters(), lr=0.01)


start_time =time.time()

for epoch in range(epochs):
    total_loss = 0
    for x_batch, y_batch in dataloader:
        x_batch, y_batch = x_batch.to(device), y_batch.to(device)

        optimizer.zero_grad()
        y_pred = model(x_batch)
        loss = criterion(y_pred, y_batch)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    if epoch % 10 == 0:
        print(f'Epoch [{epoch}/{epochs}], Loss: {total_loss / len(dataloader):.4f}')
        
        
    if (epoch % 10 == 0)  and (epoch > 0.7*epochs):
        # === Enforce Lipschitz constraint per layer ===
        with torch.no_grad():
            for layer in [model.hidden1, model.hidden2, model.output]:
                weight = layer.weight.data
                norm = torch.linalg.norm(weight, ord=2)
                scale = max(1.0, norm.item() / lam)
                print(f'The scale is calculated as [{scale}]')
                layer.weight.data = weight / scale



# Extract trained parameters for Two hidden layer model
W1 = model.hidden1.weight.detach().cpu().numpy()
b1 = model.hidden1.bias.detach().cpu().numpy()
W2 = model.hidden2.weight.detach().cpu().numpy()
b2 = model.hidden2.bias.detach().cpu().numpy()
W3 = model.output.weight.detach().cpu().numpy()
b3 = model.output.bias.detach().cpu().numpy()

# Reshape y_std to [10, 1] for correct broadcasting
y_std_np = y_std.numpy().reshape(-1, 1)

# Denormalize final layer (W3, b3)
W3_denorm = (y_std_np * W3)
b3_denorm = (y_std.numpy() * b3) + y_mean.numpy()


end_time = time.time()

training_time = end_time - start_time

# Save to MATLAB format
save_path = r"Temp_files_mid_run/trained_relu_weights_2h_norm.mat"
scipy.io.savemat(save_path, {
    'W1': W1, 'b1': b1, 'W2': W2, 'b2': b2, 'W3': W3_denorm, 'b3': b3_denorm, 'Model_training_time': training_time
})

print(f"Trained parameters saved to {save_path}")

del x, y, x_batch, y_batch, y_pred 
gc.collect()
torch.cuda.empty_cache()