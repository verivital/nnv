function Hf=hessianTensor_vanderPolEq(x,u)



 Hf{1} = interval(sparse(3,3),sparse(3,3));



 Hf{2} = interval(sparse(3,3),sparse(3,3));

Hf{2}(1,1) = -2*x(2);
Hf{2}(2,1) = -2*x(1);
Hf{2}(1,2) = -2*x(1);
