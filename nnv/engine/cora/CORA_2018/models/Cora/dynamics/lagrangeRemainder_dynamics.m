function lf=lagrangeRemainder_dynamics(x,u,dx,du)

lf=interval();

lf(1,1)=dx(3)^2/2;
lf(2,1)=0;
lf(3,1)=0;
