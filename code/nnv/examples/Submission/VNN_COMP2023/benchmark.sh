#!/bin/bash

while IFS="," read -r onnx vnnlib timeout
do
    # arguments: [1] version type, [2] benchmark category, [3] onnx file, [4] vnnlib file

    bash prepare_instance.sh 'v1' 'mnistfc' $onnx $vnnlib

    bash run_instance.sh 'v1' 'mnistfc' $onnx $vnnlib $timeout results.txt

done < examples/cifar2020_instances.csv
