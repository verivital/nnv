function dx = tanh_plant(x,u)
pl = load("plant_3rd_order_tanh.mat");
w1 = double(extractdata(pl.neuralOdeParameters.fc1.Weights));
b1 = double(extractdata(pl.neuralOdeParameters.fc1.Bias));
w2 = double(extractdata(pl.neuralOdeParameters.fc2.Weights));
b2 = double(extractdata(pl.neuralOdeParameters.fc2.Bias));
x1 = x(2,:);
x2 = x(3,:);
x4 = x(5,:);
x5 = x(6,:);
y = tanh(w1*x([3,6,7,8],:) + b1);
dx = w2*y + b2;
dx = [x1;x2;dx(1,:);x4;x5;dx(2:4,:)];
end

