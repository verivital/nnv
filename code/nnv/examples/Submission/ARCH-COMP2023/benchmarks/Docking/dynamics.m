function [dstates] = dynamics(states,actions)
    m = 12;
    n = 0.001027;

    x = states(1);
    y = states(2);
    x_dot = states(3);
    y_dot = states(4);

    Fx = actions(1);
    Fy = actions(2);

    dstates(1,1) = x_dot;
    dstates(2,1) = y_dot;
    dstates(3,1) = 2*n*y_dot + 3*(n^2)*x + Fx/m;
    dstates(4,1) = -2*n*x_dot + Fy/m;
end