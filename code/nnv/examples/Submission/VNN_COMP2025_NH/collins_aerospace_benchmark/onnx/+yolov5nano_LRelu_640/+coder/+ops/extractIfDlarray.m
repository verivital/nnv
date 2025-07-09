%#codegen
function y = extractIfDlarray(x)
if isdlarray(x)
    y = extractdata(x);
else
    y = x;
end
end
