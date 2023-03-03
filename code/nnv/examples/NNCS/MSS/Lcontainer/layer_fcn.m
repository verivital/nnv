function fcn = layer_fcn(activation_fcn)
    fcn = activation_fcn;
    if fcn == 'relu  '
        fcn = 'poslin';
    elseif fcn == 'linear'
        fcn = 'purelin';
    end
end

