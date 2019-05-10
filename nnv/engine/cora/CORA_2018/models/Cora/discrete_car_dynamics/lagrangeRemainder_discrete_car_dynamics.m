function lf=lagrangeRemainder_discrete_car_dynamics(x,u,dx,du,T)

lf=interval();

lf(1,1)=0;
lf(2,1)=(9*dx(1)^2*cos(3*x(1)))/800;
