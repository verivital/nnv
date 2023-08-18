function out = verify_spec(x)
    if class(x) ==  "Star"
        [m,M] = x.getRanges();
        x1 = max(abs(m(1)),abs(M(1)));
        x2 = max(abs(m(2)),abs(M(2)));
        x3 = max(abs(m(3)),abs(M(3)));
        x4 = max(abs(m(4)),abs(M(4)));
    elseif ismatrix(x)
        x1 = x(1); 
        x2 = x(2); 
        x3 = x(3); 
        x4 = x(4);
    else
        error('Worng input');
    end
    disp("Right side")
    disp(0.2 + 2*0.001027*sqrt(x1^2 + x2^2));
    disp("Left side")
    disp(sqrt(x3^2 + x4^2));
    if sqrt(x3^2 + x4^2) <= 0.2 + 2*0.001027*sqrt(x1^2 + x2^2)
        out = true;
    else
        out = false;
    end
end

