# Running NNV on Docker

### Requirements
 - Approx 20GB (18.2GB) of memory 
 - MATLAB license

There are two variables that can be modified to build the docker image (These are to be modified in the Dockerfile)

1. MATLAB release (default R2024b, latest supported)
`ARG MATLAB_RELEASE=R2024b`

2. License server information using the format: port@hostname
`ARG LICENSE_SERVER="27009@licesenceserver.it.vanderbilt.edu"`

Then, proceed to build the Docker image

`docker build . -t nnv`

Run Docker image interactively to use MATLAB and NNV

`docker run -it nnv`

### Acknowledgements

This work is supported in part by AFOSR, DARPA, NSF.

### Contact

For any questions related to NNV in Docker, please add them to the issues or contact [Diego Manzanas Lopez](mailto:diego.manzanas.lopez@vanderbilt.edu).

