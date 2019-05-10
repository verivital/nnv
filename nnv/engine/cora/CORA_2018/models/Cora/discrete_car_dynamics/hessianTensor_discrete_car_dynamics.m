function Hf=hessianTensor_discrete_car_dynamics(x,u,T)



 Hf{1} = interval(sparse(3,3),sparse(3,3));



 Hf{2} = interval(sparse(3,3),sparse(3,3));

Hf{2}(1,1) = (9*cos(3*x(1)))/400;
