function dx = node3(x,u)
    load("../../benchmark_dynamics/DampedOsc/ilnode(3)/5/model.mat");
    fc1 = tanh(double(w2)*x+double(b2)');
    fc2 = tanh(double(w3)*fc1+double(b3)');
    dx = double(w4)*fc2+double(b4)';
end
