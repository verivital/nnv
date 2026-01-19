# Instructions

When building the docker image, we are copying the local NNV files into the container. The result of this is that
the files are not in an editable state. So, if we make a change locally to the files, the changes will not show
up within the container. We will have to rebuild the container if we want the files to be updated.

Additionally, we need to run the `docker build` command from NNV's root so that all of the code gets copied into the
docker container correctly. We can do this with
```
cd /home/sasakis/v/tools/nnv
docker build -t nnv3.0 -f code/nnv/examples/NNV3.0/Dockerfile .
```
Then access the container with
```
docker run --gpus all -it nnv3.0
```
We need the `--gpus all` flag for the probabilistic verification which uses GPU.

# FairNNV
Make sure to run the scripts after having already run
```
cd /path/to/NNV3.0/FairNNV
```
i.e., make sure you are in the FairNNV directory before running the scripts with MATLAB headless.

After this, we can run the experiments with
```
matlab -nodisplay -r "run('run_fm26_fairnnv.m'); exit()"
```
Once the experiments run, we can copy them from the docker container onto the local filesystem by opening
a new terminal, finding the container id with `docker ps -a` and then running
```
docker cp <container_id>:/path/to/fm26_fairnnv_results /target/destination/for/fm26_fairnnv_results
```
The first path is the path to the results directory within the container and the second will be the path to the 
target destination of the directory on the local file system.

# Probabilistic Verification
After building the container and running it like in the earlier instructions, we can run the following to run the probabilistic
verification:
```
cd /path/to/NNV3.0/ProbVer
matlab -nodisplay -r "run('run_probver.m'); exit()"
```

# ModelStar

All of the code for running the examples for ModelStar are contained elsewhere within the tool (see `nnv/code/nnv/examples/Tutorial/NN/MNIST/weightPerturb`).

In this repository, run the `run_expt_for_compute.m` script to run the experiments and then produce the plot shown in the paper by running the `EXPT.m` script.

# VideoStar

The VideoStar experiments run verification on video classification neural networks using the ZoomIn dataset.
The full experiments are in `nnv/code/nnv/examples/Submission/FORMALISE2025`, and this directory contains
scripts for running a subset of those experiments (ZoomIn-4f).

## Prerequisites

The VideoStar experiments require the ONNX models and data files to be placed in the FORMALISE2025 directory structure.

### Data Setup

If you have data files in `VideoStar/data/`, copy them to the FORMALISE2025 data directory:
```bash
cp -r code/nnv/examples/NNV3.0/VideoStar/data/* code/nnv/examples/Submission/FORMALISE2025/data/
```

The expected directory structure in `FORMALISE2025/data/` is:
```
data/
└── ZoomIn/
    ├── mnistvideo_zoom_in_4f_test_data_seq.npy
    └── mnistvideo_zoom_in_test_labels_seq.npy
```

The ONNX models (e.g., `zoomin_4f.onnx`) should be in `FORMALISE2025/models/`.

## Running VideoStar ZoomIn-4f (Subset)

After building the container and running it, navigate to the VideoStar directory and run:

```bash
cd /home/matlab/nnv/code/nnv/examples/NNV3.0/VideoStar
matlab -nodisplay -r "run('run_zoomin_4f.m'); exit()"
```

Alternatively, you can use the shell script:
```bash
cd /home/matlab/nnv/code/nnv/examples/NNV3.0/VideoStar
./run_videostar_zoomin4f.sh
```

Or use the Python interface:
```bash
cd /home/matlab/nnv/code/nnv/examples/NNV3.0/VideoStar
python run_zoomin_4f.py --algorithm relax --num-samples 10
```

## Configuration

The `run_zoomin_4f.m` script can be configured by editing the CONFIGURATION section:
- `config.sampleIndices`: Which samples to verify (default: 1:10)
- `config.verAlgorithm`: Either 'relax' or 'approx' (default: 'relax')
- `config.timeout`: Timeout per sample in seconds (default: 1800)

## Results

Results are saved to `/tmp/results/VideoStar/ZoomIn/4/` with CSV files for each epsilon value:
- `eps=1_255.csv`: Results with epsilon = 1/255
- `eps=2_255.csv`: Results with epsilon = 2/255
- `eps=3_255.csv`: Results with epsilon = 3/255

To copy results from the container to your local filesystem:
```bash
docker cp <container_id>:/tmp/results/VideoStar ./videostar_results
```
