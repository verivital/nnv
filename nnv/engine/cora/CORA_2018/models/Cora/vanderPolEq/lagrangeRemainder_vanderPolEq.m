function lf=lagrangeRemainder_vanderPolEq(x,u,dx,du)

lf=interval();

lf(1,1)=0;
lf(2,1)=- dx(1)*(dx(1)*x(2) + dx(2)*x(1)) - dx(1)*dx(2)*x(1);
