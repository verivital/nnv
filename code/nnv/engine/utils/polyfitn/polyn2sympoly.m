function sp = polyn2sympoly(polyn)
% polyn2sympoly: convert a regression polynomial from polyfitn into a sympoly
% usage: sp = polyn2sympoly(polyn)
%
% arguments: (input)
%  polyn - a structure as returned from polyfitn
%
% arguments: (output)
%  sp - A sympoly object
%
% After conversion into a sympoly, any symbolic operations are
% now possible on this form.
% 
% Requirement: The sympoly toolbox, as found on Matlab Central.
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9577&objectType=FILE
%
% See also: polyvaln, polyfit, polyval, sympoly
%
% Author: John D'Errico
% Release: 1.0
% Release date: 2/19/06

if exist('sympoly','file')~=2
  error 'Please download the sympoly toolbox from Matlab Central'
end

% Copy over the fields of polyn into sp
sp.Var = polyn.VarNames;
sp.Exponent = polyn.ModelTerms;
sp.Coefficient = polyn.Coefficients(:);

% Was the list of variable names empty?
% If so, then generate a list of names of my own.
if isempty(sp.Var)
  p = size(polyn.ModelTerms,2);
  varlist={};
  for i = 1:p
    varlist{i} = ['X',num2str(i)];
  end
  sp.Var = varlist;
end

% turn the struct into a sympoly
sp = sympoly(sp);


