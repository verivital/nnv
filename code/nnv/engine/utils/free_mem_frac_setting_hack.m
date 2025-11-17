function frac = free_mem_frac_setting_hack()    % if free memory fraction is less than this, do not launch any more linear programs
% The following line MUST BE line no. 3! It will be changed by the set_free_mem_frac_for_LPs() function defined in NN.m
frac = 0.1;
end
