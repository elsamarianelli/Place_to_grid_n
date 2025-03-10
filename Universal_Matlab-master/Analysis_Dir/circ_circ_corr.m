function [circAssoc] = circ_circ_corr(theta,phi,k)
% Determine correlation between 2 circular random variables (writen by AJ)
% NB this is a NOT a rank method - cb_cc_correlation is a rank method.
%
%   CIRCULAR-CIRCULAR ASSOCIATION: T-LINEAR ASSOCIATION
%   theta and phi: circular random variables. Must be in radians
%   k: specify no. of randomisation or bootstrap samples to use in
%   estimating CIs and p-vals.
%
%   Method followed is the same as in $6.3.3 Fisher (1993), Statistical
%   Analysis of Circular Data, Cambridge University Press, ISBN: 0 521 56890 0
%
% Example data from Fisher 1992: $6.3.3, p153:
%
%   theta = [356, 97, 211, 232, 343, 292, 157, 302, 335, 302, 324, 85,...
%   324, 340, 157, 238, 254, 146, 232, 122, 329]*(pi/180);
%   phi = [119, 162, 221, 259, 270, 29, 97, 292, 40, 313, 94, 45, 47,...
%   108, 221, 270, 119, 248, 270, 45, 23]*(pi/180);
%
% NB I (CB) checked this code agains these values - it returns the correct
% Rho (0.191)
%
%   Try:
%   [circAssoc] = circ_circ_corr(theta,phi,1000)
%
% RETURNS
% circAssoc.rho     A correlation value between -1 and 1 roughtly analogous
%                   to a normal linear correlation value
%
% circAssoc.p       P value for rho
%
% circAssoc.jackknife   Conf intervals based on jackknife procedure
%
% circAssoc.bootstrap   Conf intervals based on bootstrap procedure

% Column vector inputs
theta = theta(:);
phi = phi(:);

% Test for equal sample sizes
if length(theta) ~= length(phi)
    error('length of inputs must be the same');
end

% Estimate correlation rho
[circAssoc.rho,n] = ccc(theta,phi);

% Most conservative and intuitive is to perform a randomisation test.
ind = zeros(n,k);
for i=1:k
    ind(:,i) = randperm(n)';
end
rhosim = sort(ccc(repmat(theta,1,k),phi(ind)));
[m b] = min(abs(rhosim-circAssoc.rho));
if b>k/2
    circAssoc.p = 1 - (b-1)/k;
else
    circAssoc.p = (b+1)/k;
end
clear ind rhosim

% Obtain jacknife estimates of rho and its confidence intervals
rhojack = jackknife(@ccc,theta,phi);
rhojack = n*circAssoc.rho - (n - 1)*rhojack;
circAssoc.jackknife.rho = mean(rhojack);
rhojackstd = std(rhojack);
circAssoc.jackknife.ci = [circAssoc.jackknife.rho - (1/sqrt(n))*rhojackstd*1.96,...
    circAssoc.jackknife.rho + (1/sqrt(n))*rhojackstd*1.96];

% Obtain bootstrap estimates of rho and its confidence intervals
circAssoc.bootstrap.rho = mean(bootstrp(k,@ccc,theta,phi));
circAssoc.bootstrap.ci = bootci(k,{@ccc,theta,phi})';

% --- SUB FUNCTION
function [rho,n] = ccc(theta,phi)
% Subfunction to calculate correlation between two random circular variables
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


%% These lines calculate the p-value for the correlation using conventional
%% procedures. Less computationally heavy but also less
%% intuitive/conservative than the randomisation test currently employed.
%% I've kept the code here in case it's needed for any reason in future.

% NB need to provide variable 'model' to the following lines; chooses
% hypothesis type.

%%% Might not need the following - why not just stick to randomisation
%%% test?
%
% % Test either distribution for deviation from uniformity
% [ptheta] = circ_rtest(theta);
% [pphi] = circ_rtest(phi);
% 
% % If not significantly uniform then use simple method to obtain p-value
% % otherwise use full method
% if ptheta > 0.05 || pphi > 0.05
%     if model == 0 || isempty(model)
%         circAssoc.p = exp(-abs(n*circAssoc.rho));
%     elseif model == 1
%         circAssoc.p = 0.5*exp(-n*circAssoc.rho);
%     elseif model == -1
%         circAssoc.p = 0.5*exp(n*circAssoc.rho);
%     end
% else
%     thetabar = circ_mean(theta);
%     phibar = circ_mean(phi);
%     rtheta = circ_r(theta);
%     rphi = circ_r(phi);
%     
%     alpha2theta = (1/n)*sum(cos(2*(theta - thetabar)));
%     alpha2phi = (1/n)*sum(cos(2*(phi - phibar)));
%     
%     beta2theta = (1/n)*sum(sin(2*(theta - thetabar)));
%     beta2phi = (1/n)*sum(sin(2*(phi - phibar)));
%     
%     Utheta = (1 - alpha2theta^2 - beta2theta^2)/2;
%     Uphi = (1 - alpha2phi^2 - beta2phi^2)/2;
%     Vtheta = (1 - alpha2theta)*rtheta^2;
%     Vphi = (1 - alpha2phi)*rphi^2;
%     
%     Z = (sqrt(n)*Utheta*Uphi*circAssoc.rho)/sqrt(Vtheta*Vphi);
%     
%     % Is this bit is wrong - can I be bothered to figure it out?
%     if model == 0 || isempty(model)
%         circAssoc.p = 1-normcdf(abs(Z));
%     elseif model == 1
%         circAssoc.p = 1-normcdf(Z);
%     elseif model == -1
%         circAssoc.p = normcdf(Z);
%     end
% end
