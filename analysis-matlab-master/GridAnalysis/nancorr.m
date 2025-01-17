function [ r,p ] = nancorr( var1, var2, corrType )
%NANCORR For a pair of ND arrays, remove any points with nan or inf in one
%or both arrays and then use Matlab's corr on the remaining pairs of points.

var1=var1(:);
var2=var2(:);

validPair = ~any([isnan([var1 var2]) isinf([var1, var2])] , 2);

if sum(validPair) == 0 || ~isreal(var1) || ~isreal(var2)
    r = NaN;
    p = NaN;
else
    [r,p] = corr(var1(validPair), var2(validPair), 'type', corrType);
end

end

