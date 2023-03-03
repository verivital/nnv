function isSafe = ACASXuVerifSafetyness(init_set)

% Checks if we are in a Safe scenario.
%
% If we had to stop the simulation that means there was a violation of the
% ownship's safe area (500 ft around itself) and therefore we have not
% reached a safe scenario.
%
% INPUTS
%
% init_set: state of the physical component of the system in Star form.
%
% OUTPUTS
%
% isSafe: boolean that returns True if we are in a safe state.

isSafe = true;

[x_own_lb,x_own_ub] = init_set(1).getRange(1);
[y_own_lb,y_own_ub] = init_set(1).getRange(2);
[x_int_lb,x_int_ub] = init_set(1).getRange(4);
[y_int_lb,y_int_ub] = init_set(1).getRange(5);

x_own = [x_own_lb,x_own_ub];
y_own = [y_own_lb,y_own_ub];
x_int = [x_int_lb,x_int_ub];
y_int = [y_int_lb,y_int_ub];

% calculate rho
delta_x = [x_int(1) - x_own(2), x_int(2) - x_own(1)];
delta_y = [y_int(1) - y_own(2), y_int(2) - y_own(1)];
delta_x_squared = [max(0.0, min(delta_x(1) * delta_x(1), min(delta_x(2) * delta_x(2), delta_x(1) * delta_x(2)))), ...
    max(delta_x(1) * delta_x(1), delta_x(2) * delta_x(2))];
delta_y_squared = [max(0.0, min(delta_y(1) * delta_y(1), min(delta_y(2) * delta_y(2), delta_y(1) * delta_y(2)))), ...
    max(delta_y(1) * delta_y(1), delta_y(2) * delta_y(2))];
sum_delta_squared = [delta_x_squared(1) + delta_y_squared(1), delta_x_squared(2) + delta_y_squared(2)];
rho = [sqrt(sum_delta_squared(1)),sqrt(sum_delta_squared(2))];

if rho(1) <= 500
    
    isSafe = false;
    
end
end