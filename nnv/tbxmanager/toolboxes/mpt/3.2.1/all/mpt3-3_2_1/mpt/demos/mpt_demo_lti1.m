function mpt_demo_lti1
%
% mpt_demo_lti1
%

% This demo illustrates basic usage of LTISystem. We want to model the
% discrete-time LTI system
%    x(k+1) = Ax(k) + Bu(k) 
%      y(k) = Cx(k) + Du(k)

A = [1 1; 0 1];
B = [1; 0.5];
C = [1 0];
D = 0;
lti = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

% Set the initial state of the system
lti.initialize([1; 1.5]);

% Ask for the states
disp('Initial state:');
x = lti.getStates()

% Update the system's state using some control action
disp('State update using u=0.5:');
u = 0.5;
lti.update(u); % this updates the internal state
x = lti.getStates() % ask for the updated states


% Value of the state update can also be directly obtained from update().
disp('Next state update using u=-0.6:');
u = -0.6;
next_x = lti.update(u)

disp('System output using the last known state:');
y = lti.output()

end
