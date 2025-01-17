function [new_M, R] = SR_Update(r,phi,n_phi,M,R,alpha)
%SR_UPDATE Summary of this function goes here
%   Detailed explanation goes here
gamma = 0.995; %normal is 0.995
r_alpha = 1e-6;

new_M = M + alpha * (phi' + gamma * M * n_phi' - M * phi') * phi;

%R = R + alpha * n_phi' .* (r - R);
R = R + r_alpha * M * phi' .* (r + gamma * dot(M * n_phi', R) - dot(M * phi', R));

%  [V,D] = eig(new_M);
%  V = real(V);
%  V(V<0) = 0;
%  new_M = V*D*inv(V);
end

