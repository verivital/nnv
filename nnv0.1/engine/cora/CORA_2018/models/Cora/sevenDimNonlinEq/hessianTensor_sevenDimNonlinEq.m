function Hf=hessianTensor_sevenDimNonlinEq(x,u)



 Hf{1} = sparse(8,8);



 Hf{2} = sparse(8,8);



 Hf{3} = sparse(8,8);

Hf{3}(3,2) = -4/5;
Hf{3}(2,3) = -4/5;


 Hf{4} = sparse(8,8);

Hf{4}(4,3) = -13/10;
Hf{4}(3,4) = -13/10;


 Hf{5} = sparse(8,8);

Hf{5}(5,4) = -1;
Hf{5}(4,5) = -1;


 Hf{6} = sparse(8,8);



 Hf{7} = sparse(8,8);

Hf{7}(7,2) = -3/2;
Hf{7}(2,7) = -3/2;
