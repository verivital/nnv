function dx = node3(x,u)
%     load("C:\Users\diego\Documents\GitHub\Python\sonode\experiments\damped_oscillators\ilnode(3)\5\model.mat");
    load("/home/manzand/Documents/Python/sonode/experiments/damped_oscillators/ilnode(3)/5/model.mat");
    fc1 = tanh(double(w2)*x+double(b2)');
    fc2 = tanh(double(w3)*fc1+double(b3)');
    dx = double(w4)*fc2+double(b4)';
end
