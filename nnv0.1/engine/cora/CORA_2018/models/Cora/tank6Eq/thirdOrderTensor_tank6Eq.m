function Tf = thirdOrderTensor_tank6Eq(x,u)



 Tf{1,1} = interval(sparse(7,7),sparse(7,7));

Tf{1,1}(1,1) = -897680497035489/(36028797018963968*x(1)^(5/2));


 Tf{1,2} = interval(sparse(7,7),sparse(7,7));



 Tf{1,3} = interval(sparse(7,7),sparse(7,7));



 Tf{1,4} = interval(sparse(7,7),sparse(7,7));



 Tf{1,5} = interval(sparse(7,7),sparse(7,7));



 Tf{1,6} = interval(sparse(7,7),sparse(7,7));



 Tf{1,7} = interval(sparse(7,7),sparse(7,7));



 Tf{2,1} = interval(sparse(7,7),sparse(7,7));

Tf{2,1}(1,1) = 897680497035489/(36028797018963968*x(1)^(5/2));


 Tf{2,2} = interval(sparse(7,7),sparse(7,7));

Tf{2,2}(2,2) = -897680497035489/(36028797018963968*x(2)^(5/2));


 Tf{2,3} = interval(sparse(7,7),sparse(7,7));



 Tf{2,4} = interval(sparse(7,7),sparse(7,7));



 Tf{2,5} = interval(sparse(7,7),sparse(7,7));



 Tf{2,6} = interval(sparse(7,7),sparse(7,7));



 Tf{2,7} = interval(sparse(7,7),sparse(7,7));



 Tf{3,1} = interval(sparse(7,7),sparse(7,7));



 Tf{3,2} = interval(sparse(7,7),sparse(7,7));

Tf{3,2}(2,2) = 897680497035489/(36028797018963968*x(2)^(5/2));


 Tf{3,3} = interval(sparse(7,7),sparse(7,7));

Tf{3,3}(3,3) = -897680497035489/(36028797018963968*x(3)^(5/2));


 Tf{3,4} = interval(sparse(7,7),sparse(7,7));



 Tf{3,5} = interval(sparse(7,7),sparse(7,7));



 Tf{3,6} = interval(sparse(7,7),sparse(7,7));



 Tf{3,7} = interval(sparse(7,7),sparse(7,7));



 Tf{4,1} = interval(sparse(7,7),sparse(7,7));



 Tf{4,2} = interval(sparse(7,7),sparse(7,7));



 Tf{4,3} = interval(sparse(7,7),sparse(7,7));

Tf{4,3}(3,3) = 897680497035489/(36028797018963968*x(3)^(5/2));


 Tf{4,4} = interval(sparse(7,7),sparse(7,7));

Tf{4,4}(4,4) = -897680497035489/(36028797018963968*x(4)^(5/2));


 Tf{4,5} = interval(sparse(7,7),sparse(7,7));



 Tf{4,6} = interval(sparse(7,7),sparse(7,7));



 Tf{4,7} = interval(sparse(7,7),sparse(7,7));



 Tf{5,1} = interval(sparse(7,7),sparse(7,7));



 Tf{5,2} = interval(sparse(7,7),sparse(7,7));



 Tf{5,3} = interval(sparse(7,7),sparse(7,7));



 Tf{5,4} = interval(sparse(7,7),sparse(7,7));

Tf{5,4}(4,4) = 897680497035489/(36028797018963968*x(4)^(5/2));


 Tf{5,5} = interval(sparse(7,7),sparse(7,7));

Tf{5,5}(5,5) = -897680497035489/(36028797018963968*x(5)^(5/2));


 Tf{5,6} = interval(sparse(7,7),sparse(7,7));



 Tf{5,7} = interval(sparse(7,7),sparse(7,7));



 Tf{6,1} = interval(sparse(7,7),sparse(7,7));



 Tf{6,2} = interval(sparse(7,7),sparse(7,7));



 Tf{6,3} = interval(sparse(7,7),sparse(7,7));



 Tf{6,4} = interval(sparse(7,7),sparse(7,7));



 Tf{6,5} = interval(sparse(7,7),sparse(7,7));

Tf{6,5}(5,5) = 897680497035489/(36028797018963968*x(5)^(5/2));


 Tf{6,6} = interval(sparse(7,7),sparse(7,7));

Tf{6,6}(6,6) = -897680497035489/(36028797018963968*x(6)^(5/2));


 Tf{6,7} = interval(sparse(7,7),sparse(7,7));

