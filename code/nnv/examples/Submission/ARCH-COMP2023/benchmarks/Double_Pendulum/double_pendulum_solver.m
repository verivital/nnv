syms x1 x2 th1 th2 u1 u2 g mass l c T1 T2 'real'

ex1 = 2 * x1 + x2 * cos(th2 - th1) - u2 ^ 2 * sin(th2 - th1) + ...
          - 2 * g / l * sin(th1) + (c * u1 - T1) / (mass * l ^ 2);
     
ex2 = x1 * cos(th2 - th1) + x2 + u1 ^ 2 * sin(th2 - th1) + ...
      - g / l * sin(th2) + (c * u2 - T2) / (mass * l ^ 2);
  
x1s = solve(ex1, x1);
ex2s = subs(ex2, x1, x1s);
x2ss = solve(ex2s, x2);
x1ss = subs(x1s, x2, x2ss);

x1ss_sub = subs(x1ss, [mass, g, l, c], [0.5, 1, 0.5, 0]);
%x1ss_sub = subs(x1ss, [mass, g, l, c], [1, 1, 1, 0]);
disp(x1ss_sub);
disp(simplify(x1ss_sub))


x2ss_sub = subs(x2ss, [mass, g, l, c], [0.5,1,0.5, 0]);
%x2ss_sub = subs(x2ss, [mass, g, l, c], [1,1,1, 0]);
disp(x2ss_sub);
disp(simplify(x2ss_sub))

% fprintMatPy('du1', {'th1', 'th2', 'u1', 'u2', 'g', 'mass', 'l', 'c', 'T1', 'T2'}, x1ss);
% fprintMatPy('du2', {'th1', 'th2', 'u1', 'u2', 'g', 'mass', 'l', 'c', 'T1', 'T2'}, x2ss);
