function [rho, p, rho_boot, p_shuffled, ci] = circCorrTLinear(theta,phi,k,alpha,hyp,conf)
%% Determine correlation between 2 circular random variables
%
%   CIRCULAR-CIRCULAR ASSOCIATION: T-LINEAR ASSOCIATION
%   theta and phi: circular random variables. Must be in radians
%   k: specify no. of randomisation or bootstrap samples to use in
%   estimating CIs and p-vals.
%
% Inputs:
%   phi &
%   theta: vectors of circular data (in radians) whose correlation is to be
%             measured.
%
%       k: number of permutations to use to calculate p-value from
%             randomisation and bootstrap estimation of confidence
%             intervals. Leave empty to calculate p-value analytically. (NB
%             confidence intervals will not be calculated).
%
%   alpha: hypothesis test level e.g. 0.05 (5%) 0.01 (1%)
%
%     hyp: hypothesis -1/0/1 [negatively correlated, correlated in either
%             direction, positively correlated]
%
%    conf: true or false to calculate confidence intervals via jackknife
%           (highly optimised) or bootstrap - native matlab function
%   Method followed is the same as in $6.3.3 Fisher (1993), Statistical
%   Analysis of Circular Data, Cambridge University Press, ISBN: 0 521 56890 0
%
% Example data from Fisher 1992: $6.3.3, p153:
%
%   theta = [356, 97, 211, 232, 343, 292, 157, 302, 335, 302, 324, 85,...
%   324, 340, 157, 238, 254, 146, 232, 122, 329]*(pi/180);
%   phi = [119, 162, 221, 259, 270, 29, 97, 292, 40, 313, 94, 45, 47,...
%   108, 221, 270, 119, 248, 270, 45, 23]*(pi/180);
%   k = 1000;
%   alpha = 0.05;
%   hyp = 0;
%   conf = false;
%
%   Try:
%   [rho p rho_b p_shuffled ci] = circ_circ_corr(theta,phi,k,alpha,hyp,conf)
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
%                           
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License

%% Housekeeping
theta = theta(:);
phi = phi(:);

% Test for equal sample sizes
if length(theta) ~= length(phi)
    error('length of inputs must be the same');
end

%% Estimate correlation rho
rho = ccc(theta,phi);
n = numel(theta);

%% Derive p-values
% Most conservative and intuitive is to perform a randomisation test.
if ~isempty(k)
    p_shuffled = shuffledPValue(theta,phi,rho,k,hyp);
    p = NaN;
else
    % We believe there's a problem with this. Analytic and shuffled
    % p-values do not converge. Recommend using shuffled p-values only.
    p = analyticPValue(theta, phi, rho, n, alpha, hyp);
    p_shuffled = NaN;
end

%% Estimate confidence intervals for correlation
if n>=25 && conf
    % Obtain jacknife estimates of rho and its confidence intervals
    rhojack = ccc_jackknife(theta,phi);
    rhojack = n*rho - (n - 1)*rhojack;
    rho_boot = mean(rhojack);
    rhojackstd = std(rhojack);
    ci = [rho_boot - (1/sqrt(n))*rhojackstd*norminv(alpha/2,0,1),...
        rho_boot + (1/sqrt(n))*rhojackstd*norminv(alpha/2,0,1)];
elseif n<25 && n>4 && conf && ~isempty(k)
    % Obtain bootstrap estimates of rho and its confidence intervals
    rho_boot = mean(bootstrp(k,@ccc,theta,phi));
    ci = bootci(k,{@ccc,theta,phi})';
else
    rho_boot = NaN;
    ci = [NaN NaN];
end

end

function [rho] = ccc(theta,phi)
%% Subfunction to calculate correlation between two random circular variables
n = size(theta,1);
A = sum(cos(theta).*cos(phi));
B = sum(sin(theta).*sin(phi));
C = sum(cos(theta).*sin(phi));
D = sum(sin(theta).*cos(phi));
E = sum(cos(2*theta));
F = sum(sin(2*theta));
G = sum(cos(2*phi));
H = sum(sin(2*phi));

rho = 4*(A.*B - C.*D)./sqrt((n^2 - E.^2 - F.^2).*(n^2 - G.^2 - H.^2));

end

function [rho] = ccc_jackknife(theta,phi)
%% Subfunction used to calculate jackknife estimates of correlation
n = size(theta,1) - 1;

A = cos(theta).*cos(phi);
A = sum(A) - A;

B = sin(theta).*sin(phi);
B = sum(B) - B;

C = cos(theta).*sin(phi);
C = sum(C) - C;

D = sin(theta).*cos(phi);
D = sum(D) - D;


E = cos(2*theta);
E = sum(E) - E;

F = sin(2*theta);
F = sum(F) - F;

G = cos(2*phi);
G = sum(G) - G;

H = sin(2*phi);
H = sum(H) - H;

rho = 4*(A.*B - C.*D)./sqrt((n^2 - E.^2 - F.^2).*(n^2 - G.^2 - H.^2));

end

function [p_shuffled] = shuffledPValue(theta,phi,rho,k,hyp)
%% Subfunction to calculate shuffled p-values for correlation
n = size(theta,1);
E = sum(cos(2*theta));
F = sum(sin(2*theta));
G = sum(cos(2*phi));
H = sum(sin(2*phi));

cp = cos(phi);
sp = sin(phi);

ind = uint8(randperm_multi(n,k));

thetaPerms = theta(ind);
A = cp'*cos(thetaPerms);
B = sp'*sin(thetaPerms);
C = sp'*cos(thetaPerms);
D = cp'*sin(thetaPerms);

rhosim = 4*(A.*B - C.*D)./sqrt((n^2 - E.^2 - F.^2).*(n^2 - G.^2 - H.^2));

if hyp == 1
    p_shuffled = sum(rhosim>=rho)/k;
elseif hyp == -1
    p_shuffled = sum(rhosim<=rho)/k;
elseif hyp == 0
    p_shuffled = sum(abs(rhosim)>abs(rho))/k;
else
    p_shuffled = NaN;
end

end

function [ perm ] = randperm_multi( n,k )
% produces a matrix with k columns coresponding, each of which is a 
% random permutations of 1:n
%
% DM, based on Matlab's old version of randperm
%

[ignore, perm] = sort(rand(n,k),1);

end

function p = analyticPValue(theta, phi, rho, n, alpha, hyp)
%% We can also derive p-value analytically
% Test either distribution for deviation from uniformity
[ptheta] = circRTest(theta);
[pphi] = circRTest(phi);

% If either are uniform then use simple method to get p-value
if ptheta > alpha || pphi > alpha
    if hyp == 0 || isempty(hyp)
        p = exp(-abs(n*rho));
    elseif hyp == 1
        p = 0.5*exp(-n*rho);
    elseif hyp == -1
        p = 0.5*exp(n*rho);
    end
% otherwise use full method
else
    thetabar = circMean(theta);
    phibar = circMean(phi);
    rtheta = circR(theta);
    rphi = circR(phi);
    
    alpha2theta = (1/n)*sum(cos(2*(theta - thetabar)));
    alpha2phi = (1/n)*sum(cos(2*(phi - phibar)));
    
    beta2theta = (1/n)*sum(sin(2*(theta - thetabar)));
    beta2phi = (1/n)*sum(sin(2*(phi - phibar)));
    
    Utheta = (1 - alpha2theta^2 - beta2theta^2)/2;
    Uphi = (1 - alpha2phi^2 - beta2phi^2)/2;
    Vtheta = (1 - alpha2theta)*rtheta^2;
    Vphi = (1 - alpha2phi)*rphi^2;
    
    Z = (sqrt(n)*Utheta*Uphi*rho)/sqrt(Vtheta*Vphi);
    
    if hyp == 0 || isempty(hyp)
        p = (1-normcdf(abs(Z)))*2;
    elseif hyp == 1
        p = 1-normcdf(Z);
    elseif hyp == -1
        p = normcdf(Z);
    end
end
end

function pval = circRTest(alpha)
%% From Berens' Circular Stats toolbox
r =  circR(alpha);
n = numel(alpha);
% compute Rayleigh's R (equ. 27.1)
R = n*r;
% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));

end

function mu = circMean(alpha)

mu = angle(sum(exp(1i*alpha)));

end

function R = circR(alpha)

R = -abs(sum(exp(1i*alpha)))/numel(alpha);

end

