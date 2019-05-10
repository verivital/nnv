function Hf=hessianTensor_dynamics(x,u)



 Hf{1} = interval(sparse(4,4),sparse(4,4));

Hf{1}(3,3) = 1;


 Hf{2} = interval(sparse(4,4),sparse(4,4));



 Hf{3} = interval(sparse(4,4),sparse(4,4));

