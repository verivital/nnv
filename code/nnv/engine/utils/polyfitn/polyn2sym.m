function s = polyn2sym(polyn)
% polyn2sympoly: convert a regression polynomial from polyfitn into its symbolic toolbox form
% usage: sp = polyn2sym(polyn)
%
% arguments: (input)
%  polyn - a structure as returned from polyfitn
%
% arguments: (output)
%  s - A symbolic toolbox object
%
% After conversion into a symbolic toolbox form, any symbolic operations are
% now possible on the polynomial.
% 
% Requirement: The symbolic toolbox, as supplied by the MathWorks.
% http://www.mathworks.com/products/symbolic/functionlist.html
%
% See also: polyvaln, polyfit, polyval, sym
%
% Author: John D'Errico
% Release: 3.0
% Release date: 8/23/06

if exist('sym','file')~=2
  error 'Please obtain the symbolic toolbox from the MathWorks'
end

% initialize the returned argument as symbolic
s = sym(0);

% Unpack the fields of polyn for use
Varlist = polyn.VarNames;
Expon = polyn.ModelTerms;
Coef = polyn.Coefficients;
% how many terms?
nterms = length(Coef);

% Was the list of variable names empty?
% If so, then generate a list of names of my own.
nvars = size(polyn.ModelTerms,2);
if isempty(Varlist)
  Varlist={};
  for i = 1:nvars
    Varlist{i} = ['X',num2str(i)];
  end
end

% make the vars symbolic
for i = 1:nvars
  Varlist{i} = sym(Varlist{i});
end

% build the polynomial
for i = 1:nterms
  term = sym(Coef(i));
  for j = 1:nvars
    term = term*Varlist{j}^Expon(i,j);
  end
  
  % accumulate into s
  s = s+term;
end


