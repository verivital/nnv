import sys
from typing import List
import gdown
import numpy as np
import os
import json

def create_input_bounds(c: np.ndarray, z: np.ndarray, eps: float, changable_latents_idx: List) -> np.ndarray:
    """
    Creates input bounds for the given input and epsilon.

    Returns:
        A  (c.shape[0] + z.shape[0]) x 2 tensor with the lower bounds in [..., 0] and upper bounds
        in [..., 1].
    """

    bounds = np.zeros((5, 2), dtype=np.float32)
    bounds[0, 0] = np.clip((c - eps), -1, 1)
    bounds[0, 1] = np.clip((c + eps), -1, 1)
    for i in range(4):
        if i in changable_latents_idx:
            bounds[i+1, 0] = np.clip((z[i] - eps), -1, 1)
            bounds[i+1, 1] = np.clip((z[i] + eps), -1, 1)
        else:
            bounds[i+1, 0] = z[i]
            bounds[i+1, 1] = z[i]
    
    bounds = bounds.reshape(5, 2)
    return bounds


def create_output_bounds(c: np.ndarray, eps: float) -> np.ndarray:
    """
    Creates output bounds for the given input and epsilon.

    Returns:
        A  c.shape[0] x 2 tensor with the lower bounds in [..., 0] and upper bounds
        in [..., 1].
    """

    bounds = np.zeros((1, 2), dtype=np.float32)
    bounds[0, 0] = np.clip((c - eps), -1, 1)
    bounds[0, 1] = np.clip((c + eps), -1, 1)
    
    bounds = bounds.reshape(1, 2)
    return bounds

def save_vnnlib(input_bounds: np.ndarray, output_bounds: np.ndarray, spec_path: str):

    """
    Saves the property derived as vnn_lib format.
    """

    with open(spec_path, "w") as f:

        # Declare input variables.
        f.write("\n")
        for i in range(input_bounds.shape[0]):
            f.write(f"(declare-const X_{i} Real)\n")
        f.write("\n")

        # Declare output variables.
        f.write("\n")
        f.write(f"(declare-const Y_0 Real)\n")
        f.write("\n")

        # Define input constraints.
        f.write(f"; Input constraints:\n")
        for i in range(input_bounds.shape[0]):
            f.write(f"(assert (<= X_{i} {input_bounds[i, 1]}))\n")
            f.write(f"(assert (>= X_{i} {input_bounds[i, 0]}))\n")
            f.write("\n")
        f.write("\n")

        # Define output constraints.
        f.write(f"; Output constraints:\n")
        f.write(f"(assert (or")
        f.write(f"\n")
        f.write(f"    (and (>= Y_0 {output_bounds[0, 1]}))\n")
        f.write(f"    (and (<= Y_0 {output_bounds[0, 0]}))\n")
        f.write(f"))")
        f.write("\n")

if __name__ == '__main__':
              
    try:
        np.random.seed(seed=int(sys.argv[1]))
    except (IndexError, ValueError):
        raise ValueError("Expected seed (int) to be given as command line argument")
    
    

    # Download google drive folder that contains onnx files
    onnx_folder = "onnx"
    if not os.path.exists(onnx_folder):
        url = "https://drive.google.com/drive/folders/1bh5nWgS-s5QPa2LOjQknb7Lki93eja9B"
        gdown.download_folder(url, quiet=True, use_cookies=False)
    files = os.listdir(onnx_folder)
    assert len(files) == 8, "Expected 8 onnx files in the onnx_files directory"

    # create a folder for the vnnlib files
    vnnlib_folder = "vnnlib"
    if not os.path.exists(vnnlib_folder):
        os.mkdir(vnnlib_folder)

    csv_path = "instances.csv"
    f = open(csv_path, "w")

    json_file = open("parameters.json", "r")
    parameters = json.load(json_file)
    for parameter in parameters['models']:
        file = parameter['name']
        assert os.path.exists(os.path.join(onnx_folder, file)), f"Expected {file} to be in the onnx_files directory"

        num_instances = parameter['num_instances']
        input_epsilons = parameter['input_epsilons']
        num_latent_sets = parameter['num_latent_sets']
        output_epsilons = parameter['output_epsilons']
        assert len(input_epsilons) == len(output_epsilons), "Expected input_epsilons and output_epsilons to have the same length"
        timeout = parameter['timeout']

        for idx in range(num_instances):
            epsilon_idx = np.random.choice(range(len(input_epsilons)))
            input_epsilon = input_epsilons[epsilon_idx]
            output_epsilon = output_epsilons[epsilon_idx]
            latent_sets_idx = \
                np.random.choice(range(4), num_latent_sets[epsilon_idx], replace=False)
            
            c = np.random.uniform(low=0.2, high=0.7, size=(1, ))
            z = np.random.randn(4, )
            z = np.clip(z, -1+input_epsilon, 1-input_epsilon)

            input_bounds = create_input_bounds(c, z, input_epsilon, latent_sets_idx)
            output_bounds = create_output_bounds(c, output_epsilon)

            vnnlib_file_name = f"{vnnlib_folder}/{file[:-5]}_prop_{idx}_input_eps_{input_epsilon:.3f}_output_eps_{output_epsilon:.3f}.vnnlib"

            save_vnnlib(input_bounds, output_bounds, vnnlib_file_name)

            net = f"{onnx_folder}/{file}"
            f.write(f"{net},{vnnlib_file_name},{timeout}\n")
        
    f.close()
