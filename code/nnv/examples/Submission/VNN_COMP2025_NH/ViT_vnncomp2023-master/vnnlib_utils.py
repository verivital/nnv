import torch


def create_input_bounds(img: torch.Tensor, eps: float,
                        mean: torch.Tensor, std: torch.Tensor) -> torch.Tensor:
    mean = mean.view(-1, 1, 1)
    std = std.view(-1, 1, 1)
    bounds = torch.zeros((*img.shape, 2), dtype=torch.float32)
    bounds[..., 0] = (torch.clamp((img - eps), 0, 1) - mean) / std
    bounds[..., 1] = (torch.clamp((img + eps), 0, 1) - mean) / std
    return bounds.view(-1, 2)


def save_vnnlib(input_bounds: torch.Tensor, label: int, spec_path: str, total_output_class: int):
    with open(spec_path, "w") as f:

        f.write(f"; Property with label: {label}.\n")

        # Declare input variables.
        f.write("\n")
        for i in range(input_bounds.shape[0]):
            f.write(f"(declare-const X_{i} Real)\n")
        f.write("\n")

        # Declare output variables.
        f.write("\n")
        for i in range(total_output_class):
            f.write(f"(declare-const Y_{i} Real)\n")
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

        # disjunction version:
        f.write("(assert (or\n")
        for i in range(total_output_class):
            if i != label:
                f.write(f"    (and (>= Y_{i} Y_{label}))\n")
        f.write("))\n")
