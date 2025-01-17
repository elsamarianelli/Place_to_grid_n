function [t] = vonMisesRand(siz,m,k)
%% Notes
% siz - size of output matrix of von Mises random numbers
% m - size siz vector or scalar containing value(s) of the mean parameter
%     of von Mises distribution from which to draw values.
% k - size siz matrix or scalar containing value(s) of the concentration
%     parameter of von Mises distribution from which to draw values.
% 
% Example:
% x = -1:0.0001:1;
% k = exp(log(4) + log(2)*x); % vector of varying dispersion or choose
%                             % scalar value.
% m0 = pi;                    % Set model parameters for varying mean m.
% b = -pi;                    % Set model parameters for varying mean m.
% m = m0 + 2*atan(b*x);       % vector of varying mean or choose scalar
%                             % fixed value.
%
% siz = size(m) or size(k) or user defined for scalar m,k
% t = vonMisesRand(siz,m,k);
%
% Simulated values may be used to test circ_regress function for estimating
% parameters of circular-linear regression curve. Parameters, b, m0 etc are
% analogous to gradient and intercept and are the parameters of a linear
% model of how the mean parameter of a von Mises distribution varies with a
% single linear covariate x of circular response variable t.
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk>
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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

if numel(k)==1
    k = k*ones(siz);
end
if numel(m)==1
    m = m*(ones(siz));
end
if prod(siz)~=numel(m) || prod(siz)~=numel(k)
    error('mismatch in dimensions of m and k')
end

t = zeros(siz);

tau = 1 + sqrt((1 + 4*k.^2));
rho = (tau - sqrt((2*tau)))./(2.*k);
r = (1 + rho.^2)./(2*rho);

accepted = zeros(siz);
while sum(accepted(:)) ~= numel(accepted)
    u1 = rand(siz);
    u2 = rand(siz);
    u3 = rand(siz);
    
    z = cos(pi*u1);
    f = (1 + r.*z)./(r + z);
    c = k.*(r - f);
    
    accept = (log(c./u2) + 1 - c) > 0;
    ind = accept & ~accepted;
    t(ind) = sign(u3(ind) - 0.5).*acos(f(ind));
    
    accepted = accept | accepted;
end
    
t = fixAngle(t + m);

end