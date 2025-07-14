import os
import argparse
import csv
import torch

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('seed', default=0, type=int)
    return parser.parse_args()

def set_seed(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

def generate_text_specs(C_mat, rhs_mat):
    nonzero_cumsum = (C_mat != 0).int().cumsum(dim=-1)
    first_var_idx = (nonzero_cumsum == 1).max(dim=-1)[1]
    second_var_idx = (nonzero_cumsum == 2).max(dim=-1)[1]
    var_var_mask = (second_var_idx != 0)

    rhs_mat = torch.where(var_var_mask, torch.zeros_like(rhs_mat), rhs_mat)

    first_var_sign = C_mat.gather(2, first_var_idx.unsqueeze(-1)).squeeze(-1).sign()
    smaller_eq_mask = (first_var_sign == 1)
    C_mat = C_mat * first_var_sign.unsqueeze(-1)
    rhs_mat = rhs_mat * first_var_sign

    second_var_pos_mask = C_mat.gather(2, second_var_idx.unsqueeze(-1)).squeeze(-1) > 0

    spec_gen_list = []
    for i in range(C_mat.shape[0]):
        spec_gen = []
        for j in range(C_mat.shape[1]):
            if var_var_mask[i, j]:
                spec_text = f"({'<' if smaller_eq_mask[i, j] else '>'}= Y_{first_var_idx[i, j]} {'-' if second_var_pos_mask[i, j] else ''}Y_{second_var_idx[i, j]})"
            else:
                spec_text = f"({'<' if smaller_eq_mask[i, j] else '>'}= Y_{first_var_idx[i, j]} {rhs_mat[i, j]})"
            spec_gen.append(spec_text)
        spec_gen_list.append(spec_gen)

    return spec_gen_list

def generate_vnnlib_single_input(vnnlib_path, data_min, data_max, C_mat, rhs_mat, flat_cond_mat):
    num_distinct_X = data_min.shape[0]
    assert num_distinct_X == 1, "Assume only one input range"
    data_min = data_min.view(-1)
    data_max = data_max.view(-1)
    input_dim = data_min.shape[0]
    output_dim = C_mat.shape[-1]
    num_or = flat_cond_mat.shape[0]

    spec_list = generate_text_specs(C_mat, rhs_mat)[0]

    with open(vnnlib_path, "w") as f:
        # Declare input variables.
        for i in range(input_dim):
            f.write(f"(declare-const X_{i} Real)\n")
        f.write("\n")

        # Declare output variables.
        f.write("\n")
        for i in range(output_dim):
            f.write(f"(declare-const Y_{i} Real)\n")
        f.write("\n")

        # Define input constraints.
        f.write(f"; Input constraints:\n")
        for i in range(input_dim):
            f.write(f"(assert (<= X_{i} {data_max[i]}))\n")
            f.write(f"(assert (>= X_{i} {data_min[i]}))\n")
            f.write("\n")
        f.write("\n")

        # Define output constraints.
        f.write(f"; Output constraints:\n")

        # disjunction version:
        f.write("(assert (or\n")
        spec_idx = 0
        for i in range(num_or):
            f.write("    (and")
            for j in range(flat_cond_mat[i]):
                f.write(f" {spec_list[spec_idx]}")
                spec_idx += 1
            f.write(")\n")
        f.write("))\n")

    print(f"Generated VNNLIB file: {vnnlib_path}")

def main():
    args = get_args()
    torch.set_default_dtype(torch.float64)
    device = torch.device("cpu")
    print(f"Using device: {device}")

    bench_dir = "."
    onnx_dir_name = "onnx"
    onnx_dir = os.path.join(bench_dir, onnx_dir_name)
    assert os.path.exists(onnx_dir)

    vnnlib_dir_name = "vnnlib"
    vnnlib_dir = os.path.join(bench_dir, vnnlib_dir_name)
    os.makedirs(vnnlib_dir, exist_ok=True)

    matadata_dir_name = "matadata"
    matadata_dir = os.path.join(bench_dir, matadata_dir_name)
    assert os.path.exists(matadata_dir)

    num_vnnlib_max = 50
    timeout = 150

    onnx_files = [f for f in os.listdir(onnx_dir) if f.endswith('.onnx')]

    csv_items = []
    for onnx_file in onnx_files:
        model_name = onnx_file.split('.')[0]
        print(f"Processing model: {model_name}")

        csv_model_path = os.path.join(onnx_dir_name, onnx_file)
        matadata_path = os.path.join(matadata_dir, f"{model_name}.pt")
        matadata = torch.load(matadata_path, map_location=device)
        data_min = matadata['data_min'].to(device)
        data_max = matadata['data_max'].to(device)
        spec_mat = matadata['spec_mat'].to(device)
        rhs_mat = matadata['rhs_mat'].to(device)
        num_spec_gen = rhs_mat.shape[-1]

        set_seed(args.seed)

        num_vnnlib = min(num_vnnlib_max, rhs_mat.shape[0])
        spec_indices = torch.randperm(rhs_mat.shape[0], device=device)[:num_vnnlib].tolist()

        for i, idx in enumerate(spec_indices):
            vnnlib_name = f"{model_name}_{i}.vnnlib"
            csv_vnnlib_path = os.path.join(vnnlib_dir_name, vnnlib_name)
            vnnlib_path = os.path.join(vnnlib_dir, vnnlib_name)
            generate_vnnlib_single_input(
                vnnlib_path, data_min[idx:idx+1], data_max[idx:idx+1], spec_mat[idx:idx+1], rhs_mat[idx:idx+1], torch.ones((1,), dtype=torch.int, device=device) * num_spec_gen
            )
            csv_items.append((csv_model_path, csv_vnnlib_path, timeout))

    csv_path = os.path.join(bench_dir, "instances.csv")
    with open(csv_path, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(csv_items)

    print(f"Generated CSV file: {csv_path}")

if __name__ == "__main__":
    main()
