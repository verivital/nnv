% construct a Average Pooling 2D Layer object

L1 = MaxPooling2DLayer('test_max_pooling_2d_layer', [2 2], [1 1], [0 0 0 0]);
L2 = MaxPooling2DLayer();
L3 = MaxPooling2DLayer([3 3], [1 1], [0 0 0 0]);