function dx = node4(x,u)
    load("../../benchmark_dynamics/DampedOsc/ilnode(4)/5/model.mat");
    fc1 = tanh(double(Wb{3})*x+double(Wb{4}'));
    fc2 = tanh(double(Wb{5})*fc1+double(Wb{6}'));
    dx = double(Wb{7})*fc2+double(Wb{8})';
end
