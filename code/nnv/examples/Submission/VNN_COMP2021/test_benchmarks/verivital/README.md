### VNN Comp 2021: Proposed Benchmarks

The benchmarks consist of two MNIST classifiers, one with a maxpooling layer and the other with an average pooling layer. 

For both the networks we suggest 20(randomly selected) images:
  1. for the average pooling network with a perturbation radii of 0.02 and 0.04 and a timeout of 5 minutes.
  2. for the max pooling network with a perturbation radii of 0.004 and a timeout of 7 minutes.

The network expects the input image to be a Tensor of size (Nx1x28x28)[i.e NCHW] and normalized between [0,1].

#### Generating vnnlib specifications:

Run: 
```
python generate_properties.py --seed <seed>
```
Here 'seed' is an int variable and the default is set to 0. Based on the seed value a specific set of images will be produced.

The script generate_properties.py generates the vnnlib properties as well as the 
mnist_Conv_pool_instances.csv provides the network,specs, timeout instances. 

#### Network Performace(accuracies) on Test data:

Convnet_avgpool.onnx:  97.95 %  (Training accuracy: 98.31 %)     
Convnet_maxpool.onnx:  97.92 %  (Training accuracy: 98.75 %)

For any queries, please feel free to contact: neelanjana.pal@vanderbilt.edu
                                          or  neelanjana314@gmail.com 

Acknowledgement for getting ideas on spec generating file: Patrick Henriksen(@pat676)

