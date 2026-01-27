import torch
import numpy as np
import scipy.io
import gc
import time
import argparse

# Define the function to compute covariance
@torch.jit.script
def compute_f(A, X_batch):
    """
    Compute the function f(A) = (1 / (N * norm(A)^2)) * sum_i (A' * (X_i - mean(X)))^2
    """
    norm_A2 = torch.norm(A, p=2)
    mean_X = torch.mean(X_batch, dim=1, keepdim=True)  
    N = X_batch.shape[1]
    dX = X_batch - mean_X
    A_X = torch.matmul(A.T, dX)  # A' * (X - mean)
    
    # Compute covariance
    cov_batch = (1 / (norm_A2)) * torch.sqrt( torch.sum(torch.pow(A_X, 2)) / N )
    return cov_batch


def deflation(X, n, m, device , batch_size):
    A_list = []  # To store all principal directions
    round = -1
    Largest = float('inf')
    
    while True:
        round += 1
        # Initialize A for this round
        
        if round>0:
            del X_batch
            del X_proj
        gc.collect()
        
        A = torch.randn(n, 1, requires_grad=True, device=device)  # Trainable parameter (n x 1)
        
        Initial_lr = 100*n
        
        optimizer = torch.optim.SGD([A], lr=Initial_lr)
        
        Cov_prev = float('inf')
        
        iter = 0
        Cov = compute_f(A, X)
        print(f"Initial guess for covaiance is: {Cov.item()}")
        
        # Train A to find the principal direction for this round
        while True:
            optimizer.zero_grad()  # Zero gradients
            
            
            # Sample a random batch from X
            indices = torch.randint(0, m, (batch_size,), device=device)  
            X_batch = X[:, indices]  # Select batch
            
            # Compute the function value f(A) and its gradients
            cov_batch = compute_f(A, X_batch)
            
            # Perform gradient ascent to maximize cov_batch
            (-cov_batch).backward(retain_graph=True)  # Equivalent to maximizing cov_batch
            
            optimizer.step()  # Update A
            
            
            # Check convergence
            if iter % 100 == 0:
                Cov = compute_f(A, X)
                print(f"Round {round+1} - Iteration {iter} - Covariance: {Cov.item()}")
                
                # Check if convergence condition is met
                if round == 0:
                    if abs(Cov - Cov_prev) < Cov_prev * 1e-2:
                        print(f"Convergence achieved at iteration {iter}. Stopping optimization.")
                        break
                else:
                    if abs(Cov - Cov_prev) < Largest * 1e-2:
                        print(f"Convergence achieved at iteration {iter}. Stopping optimization.")
                        break
                    
                
                Cov_prev = Cov  # Update the previous covariance value
            
            iter += 1
        
        # Store the current principal direction A
        A = A / torch.norm( A , p=2 )
        A_list.append(A.clone().detach())
        
        # Remove the component of the data in the direction of A
        X_proj = torch.matmul(A.T, X) * A  # Project the data onto A
        del A
        print(X_proj.shape)  # Should be (64*84*11, 5*2000)
        X = X - X_proj  # Subtract projection from X to remove the direction
        X = X.detach()
        
        
        if round == 0:
            Largest = Cov
        
        
        if Cov<=Largest/100:
             print(f"A sufficient amount of principal directions ( {round+1} unit vectors) is collected.")
             break
         
    return A_list, X



def main(mat_file_path: str, num_files: int, N_dir: int,
         batch_size: int, height: int, width: int, n_class: int):
    
    mat_data = scipy.io.loadmat(mat_file_path + '/Direction_data.mat')

    Y = {}
    new_shape = (height * width * n_class, N_dir // num_files)
    for i in range(1, num_files + 1):
        key = f'Y{i}'
        YY = mat_data[key];
        YY = np.reshape(YY, new_shape, order='F')
        Y[key] = torch.tensor(YY, dtype=torch.float32)
        #Y[key] = torch.tensor(mat_data[key], dtype=torch.float32).reshape(new_shape)

    del mat_data
    gc.collect()

    X = torch.cat([Y[f'Y{i}'] for i in range(1, num_files + 1)], dim=1)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    X = X.to(device)
    
    del Y
    gc.collect()

    print("X shape:", X.shape)

    n = X.shape[0]
    m = X.shape[1]
    start_time = time.time()
    A_list, X_updated = deflation(X, n, m, device, batch_size)
    end_time = time.time()

    Training_time = end_time - start_time

    print("Updated X shape:", X_updated.shape)
    print("Principal directions found:")
    for i, A in enumerate(A_list):
        print(f"A{i+1}: {A}")
    
    path = mat_file_path + '/directions.mat'
    A_list_np = np.hstack([A.to('cpu').numpy().astype(np.float32) for A in A_list])
    scipy.io.savemat(path, {'Directions': A_list_np, 'Direction_Training_time': Training_time})
    print("Directions and training time saved to 'directions.mat'")

    del A_list, A_list_np, X
    gc.collect()
    torch.cuda.empty_cache()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Principal Direction Extraction via Deflation")
    parser.add_argument('--mat_file_path', type=str, required=True, help='Path to the input .mat file.')
    parser.add_argument('--num_files', type=int, required=True, help='Number of .mat file variables (Y1, Y2, ..., Yn).')
    parser.add_argument('--N_dir', type=int, required=True, help='Total number of data instances (divided across num_files).')
    parser.add_argument('--batch_size', type=int, required=True, help='Batch size used in optimization.')
    parser.add_argument('--height', type=int, required=True, help='Image height (e.g., 720).')
    parser.add_argument('--width', type=int, required=True, help='Image width (e.g., 960).')
    parser.add_argument('--n_class', type=int, required=True, help='Number of classes (e.g., 12).')
    
    args = parser.parse_args()

    main(args.mat_file_path, args.num_files, args.N_dir,
         args.batch_size, args.height, args.width, args.n_class)
