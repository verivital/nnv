function lf=lagrangeRemainder_accSysEidFast(x,u,dx,du)

lf=interval();

lf(1,1)=0;
lf(2,1)=- dx(2)*((1281*du(1))/(50*x(2)^2) - (1281*dx(2)*u(1))/(25*x(2)^3)) - (1281*du(1)*dx(2))/(50*x(2)^2);
