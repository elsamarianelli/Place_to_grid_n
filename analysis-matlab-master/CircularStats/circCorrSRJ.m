function [rho, pval, ts] = circCorrSRJ(alpha1, alpha2)

% [rho, p_analytic, p_shuffled] = circCorrSRJ(alpha1, alpha2,k,hyp)
%
% [rho pval ts] = circ_corrcc(alpha1, alpha2)
%   Circular correlation coefficient for two circular random variables.
%
%   Input:
%     alpha1 sample of angles in radians
%     alpha2 sample of angles in radians
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%   Topics in circular statistics, S.R. Jammalamadaka et al., p. 176
%
% PHB 6/7/2008
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
% Distributed under Creative Commons license
% Attribution-Noncommercial-Share Alike 3.0
% http://creativecommons.org/licenses/by-nc-sa/3.0/
%
%% Shuffled p-value + optimisation added by Daniel Manson 

if size(alpha1,2) > size(alpha1,1)
alpha1 = alpha1';
end

if size(alpha2,2) > size(alpha2,1)
alpha2 = alpha2';
end

if length(alpha1)~=length(alpha2)
  error('Input dimensions do not match.')
end

% compute mean directions
n = length(alpha1);
a = sin(alpha1(:) - circMean(alpha1));
b = sin(alpha2(:) - circMean(alpha2));

% compute correlation coeffcient from p. 176
num = a' * b;
den = sqrt(sum(a.^2)  * sum(b.^2));
rho = num / den; 

% compute pvalue analytically
l20 = sum(a.^2)/n;
l02 = sum(b.^2)/n;
l22 = a'.^2 * b.^2 / n;

ts = sqrt((n * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));


%% Uncomment for shuffled p-value
% % compute pvalue by shuffling
% num_shuffle = a' * b(randperm_multi(n,k));
% rho_shuffle = num_shuffle / den;
% 
% if hyp == 1
%     p_shuffled = sum(rho_shuffle>=rho)/k;
% elseif hyp == -1
%     p_shuffled = sum(rho_shuffle<=rho)/k;
% elseif hyp == 0
%     p_shuffled = sum(abs(rho_shuffle)>abs(rho))/k;
% else
%     p_shuffled = NaN;
% end

end

%% Uncomment for shuffled p-value
% function [ perm ] = randperm_multi( n,k )
% % produces a matrix with k columns coresponding, each of which is a 
% % random permutations of 1:n
% %
% % DM, based on Matlab's old version of randperm
% %
% 
% [ignore, perm] = sort(rand(n,k),1);
% 
% end

function mu = circMean(alpha)

mu = angle(sum(exp(1i*alpha)));

end
