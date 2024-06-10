function dx = dynamics(x,u)

dx = [
    x(3)*cos(x(4));
    x(3)*sin(x(4));
    u(1);
    u(2)
];

end

