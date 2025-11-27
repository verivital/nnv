function pout = polydern(p,diffvar)
% polydern: partial derivatives of a n-dimensional polynomial model created by polyfitn.
% usage: pout = polydern(p,diffvar)
%
% arguments: (input)
%  p - a structure containing a polynomial regression model.  See polyfitn
%      for details.
%
%  diffvar - integer, which variable to compute the derivative with respect
%      to. If diffvar is a character string, then it refers to an explicit
%      variable in the model that is found in the VarNames field. If so,
%      the variable name must be an exact match for one of the variables
%      listed in p.VarNames.
%
% arguments: (output)
% pout - a structure containing a polynomial regression model which
%           describes the partial derivative with respect to variable #i.
%
% Note:
% Polydern modifies the "ModelTerms" and "Coefficients" fields in creating
% the output polynomial. ParameterVar, parameterStd are also updated to
% have the proper shape, as well as scaling them to reflect the derivative.
% 
% The other fields (R2, RMSE, VarNames) are copied without change, even
% though the error measures are no longer meaningful for the derivative
% polynomial.
%
% This function uses polynomial structure variables created by polyfitn, 
% which is part of the PolyfitnTools package by John D'Errico:
%     http://www.mathworks.com/matlabcentral/fileexchange/10065
%
% See also: polyder, polyfitn, polyvaln
%
% Polydern author: Jason Goodman
% Error checks, ParameterStd/Var updates and allowance for the variable
% name itself instead of an index added by John D'Errico

% check the number of input args
if nargin ~= 2
  error('POLYDERN:argumentcount','Exactly two arguments are required')
end

% is p a struct, created by polyfitn?
if ~isstruct(p) || ~isfield(p,'ModelTerms') || ~isfield(p,'Coefficients')
  error('POLYDERN:invalidp','p must be a struct as created by polyfitn')
end

% check diffvar too. scalar, numeric, integer, real, in the proper range,
% or it must be character, an exact match for one of the variables listed
% in p.VarNames.
if ischar(diffvar)
  ind = ismember(p.VarNames,diffvar);
  if sum(ind) == 1
    % we have a hit
    diffvar = find(ind);
  else
    error('POLYDERN:invaliddiffvar','diffvar must be a scalar numeric real integer value, or a valid variable name')
  end
elseif ~isnumeric(diffvar) || (numel(diffvar) ~= 1) || (round(diffvar) ~= diffvar) || ~isreal(diffvar)
  % diffvar must be scalar, numeric, etc. Also verify that it is real to
  % be complete
  error('POLYDERN:invaliddiffvar','diffvar must be a scalar numeric real integer value, or a valid variable name')
elseif (diffvar < 1) || (diffvar > size(p.ModelTerms,2))
  % test that diffvar is at least 1, and is not larger than the
  % number of variables in the nmodel
  error('POLYDERN:invaliddiffvar', ...
    ['diffvar (=',num2str(diffvar), ...
    ') cannot represent the variable index to be differentiated, as it is out of range for this model'])
end

pout = p;
pout.ModelTerms = [];
pout.Coefficients = [];

remainingterms = true(1,size(p.ModelTerms,1));
jout = 1;
for jin=1:length(p.Coefficients)
    if (p.ModelTerms(jin,diffvar) ~= 0)
        pout.Coefficients(jout) = p.Coefficients(jin).*p.ModelTerms(jin,diffvar);
        pout.ModelTerms(jout,:) = p.ModelTerms(jin,:);
        pout.ModelTerms(jout,diffvar) = pout.ModelTerms(jout,diffvar)-1;
        jout = jout + 1;
    else
        % this term got dropped from the model
        remainingterms(jin) = false;
    end
end

if (jout == 1)   % Polynomial has no terms in this variable: deriv is zero
    [ncoeff, nvars] = size(p.ModelTerms);
    pout.ModelTerms(1,:) = zeros(1,nvars);
    pout.Coefficients = 0;
    
    pout.ParameterVar = 0;
    pout.ParameterStd = 0;
else
    % at least some terms remain, so update the parameter variances and
    % standard deviations. First, drop those terms that died.
    pout.ParameterVar = reshape(p.ParameterVar(remainingterms),1,[]);
    pout.ParameterStd = reshape(p.ParameterStd(remainingterms),1,[]);
    
    % scale them appropriately. Thus if std(A) is sigma, then std(k*A)
    % is k*sigma, but multiply the corresponding variances by k.^2.
    scale = p.ModelTerms(remainingterms,diffvar).';
    pout.ParameterStd = pout.ParameterStd.*scale;
    pout.ParameterVar = pout.ParameterVar.*scale.^2;
end

