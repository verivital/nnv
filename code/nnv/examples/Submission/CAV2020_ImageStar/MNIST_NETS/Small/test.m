outputSet = nnvNet.reachSet{14};

S1 = outputSet(1);
labels = correct_labels{1};
id1 = labels(1);

b = CNN.isRobust(S1, id1);

[lb1, ub1] = S1.getRanges;

classified_id = CNN.classifyOutputSet(S1);


b2 = S1.is_p1_larger_p2([1 1 2], [1 1 1])
b3 = S1.is_p1_larger_p2([1 1 3], [1 1 1])
b4 = S1.is_p1_larger_p2([1 1 4], [1 1 1])
b5 = S1.is_p1_larger_p2([1 1 5], [1 1 1])
b6 = S1.is_p1_larger_p2([1 1 6], [1 1 1])
b7 = S1.is_p1_larger_p2([1 1 7], [1 1 1])
b8 = S1.is_p1_larger_p2([1 1 8], [1 1 1])
b9 = S1.is_p1_larger_p2([1 1 9], [1 1 1])
b10 = S1.is_p1_larger_p2([1 1 10], [1 1 1])