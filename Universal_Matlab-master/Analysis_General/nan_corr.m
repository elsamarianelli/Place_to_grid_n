function [ r,p ] = nan_corr( var1, var2, corrType)
%NAN_CORR Corr that ignores nans, infs, and complex
%  Take variables, removes nans, infs and comples then passes to matlabs
%  built in corr function. Only works on vector inputs - not mat.

var1    =var1(:);                                                           %Make sure it's column
var2    =var2(:);

nonValNan   =any(isnan([var1, var2]),2);                                    %Get inds for nan, inf, complex
nonValInf   =any(isinf([var1, var2]),2);
nonValComplex =any(~isreal([var1, var2]),2);

validPair   =~nonValNan & ~nonValInf & ~nonValComplex;                      %Only process valid values

% If insufficent values are left then just return nan values
if sum(validPair) == 0
    r       =NaN;
    p       =NaN;
else
    [r,p]   =corr(var1(validPair), var2(validPair), 'type', corrType);
end

end

