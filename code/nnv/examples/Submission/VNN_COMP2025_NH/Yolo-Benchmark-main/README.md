# Yolo-Benchmark
This is a benchmark of Tiny YOLOv2 for VNNCOMP 2023.

## Run
`python generate_instances.py <seed>`

## Model
Our model is a modified version of YOLOv2. For simplicity of verification, the backbone of the original YOLOv2 model is replaced with a much smaller network. The network takes in an image (3 x 52 x 52) and outputs a tensor of a shape 125 x 13 x 13, which contains the confidence scores, the classification result, and the positions of 5 bounding boxes for each of the 13 x 13 grid cells.

## Specifications
Robustness properties are defined on the inputs and raw output tensors of the backbone network. For each instance, the property is verified if there is at least one bounding box successfully predicting the same object (with a confidence score greater than a preset threshold) given any bounded random perturbations on the input image.
