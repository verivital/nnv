% construct a Average Pooling 2D Layer object

L1 = AveragePooling2DLayer('test_average_pooling_2d_layer', [2 2], [1 1], [0 0 0 0]);
L2 = AveragePooling2DLayer();
L3 = AveragePooling2DLayer([3 3], [1 1], [0 0 0 0]);