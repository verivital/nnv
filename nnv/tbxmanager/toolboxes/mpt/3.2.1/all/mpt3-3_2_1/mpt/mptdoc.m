function mptdoc
%
% display documentation for MPT
%

% check the version 
v = str2double(strtok(version,'.'));

% locate mptdoc
f=which('mpt.html');

if isempty(f)
    error('The documentation for MPT is not installed on the Matlab path.');
end
if v<8    
    web(f,'-helpbrowser');
else
    % new version
    doc -classic;
end


end