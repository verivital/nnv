function lf=lagrangeRemainder_sevenDimNonlinEq(x,u,dx,du)

lf=interval();

lf(1,1)=0;
lf(2,1)=0;
lf(3,1)=-(4*dx(2)*dx(3))/5;
lf(4,1)=-(13*dx(3)*dx(4))/10;
lf(5,1)=-dx(4)*dx(5);
lf(6,1)=0;
lf(7,1)=-(3*dx(2)*dx(7))/2;
