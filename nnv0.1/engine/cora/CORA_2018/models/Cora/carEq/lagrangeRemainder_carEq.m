function lf=lagrangeRemainder_carEq(x,u,dx,du)

lf=interval();

lf(1,1)=0;
lf(2,1)=-dx(2)^2/10;
