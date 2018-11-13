
n = length(R); 

R_bounded = [];
R_unbounded = [];

for i=1:n
    fprintf('\nChecking R(%d)', i);
    if ~R(i).isBounded
        fprintf('\nR(%d) is not bounded', i);
        R_unbounded = [R_unbounded R(i)];
    else
        R_bounded = [R_bounded R(i)];
    end
end

B = [];
m = length(R_bounded);
for i=1:m
    B = [B R_bounded(i).outerApprox];
end



