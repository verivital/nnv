#!/bin/sh


chmod +x prepare_instance.sh
chmod +x run_instance.sh


# verivital
./prepare_instance.sh  'v1' 'mnist' '../test_benchmarks/verivital/Convnet_avgpool.onnx' '../test_benchmarks/verivital/specs/avgpool_specs/prop_0_0.02.vnnlib'

./run_instance.sh  'v1' 'mnist' '../test_benchmarks/verivital/Convnet_avgpool.onnx' '../test_benchmarks/verivital/specs/avgpool_specs/prop_0_0.02.vnnlib' '../results/verivital/' 300


#oval21
./prepare_instance.sh  'v1' 'cifar' '../test_benchmarks/oval21/nets/cifar_base_kw.onnx' '../test_benchmarks/oval21/vnnlib/cifar_base_kw-img4537-eps0.012679738562091505.vnnlib'

./run_instance.sh  'v1' 'cifar' '../test_benchmarks/oval21/nets/cifar_base_kw.onnx' '../test_benchmarks/oval21/vnnlib/cifar_base_kw-img4537-eps0.012679738562091505.vnnlib' '../results/oval21/' 720

#cifar2020
./prepare_instance.sh  'v1' 'cifar' '../test_benchmarks/cifar2020/nets/cifar10_2_255_simplified.onnx' '../test_benchmarks/cifar2020/specs/cifar10/cifar10_spec_idx_0_eps_0.00784_n1.vnnlib'

./run_instance.sh  'v1' 'cifar' '../test_benchmarks/cifar2020/nets/cifar10_2_255_simplified.onnx' '../test_benchmarks/cifar2020/specs/cifar10/cifar10_spec_idx_0_eps_0.00784_n1.vnnlib' '../results/cifar2020/' 100

#eran
./prepare_instance.sh  'v1' 'mnist' '../test_benchmarks/eran/nets/mnist_relu_9_200.onnx' '../test_benchmarks/eran/specs/mnist/mnist_spec_idx_1989_eps_0.01200.vnnlib'

./run_instance.sh  'v1' 'mnist' '../test_benchmarks/eran/nets/mnist_relu_9_200.onnx' '../test_benchmarks/eran/specs/mnist/mnist_spec_idx_1989_eps_0.01200.vnnlib' '../results/eran/' 300
