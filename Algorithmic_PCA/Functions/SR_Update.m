function [new_M, R] = SR_Update(r, phi, n_phi, M, R, alpha)
% SR_UPDATE Updates the Successor Representation (SR) matrix and reward estimates.
%
% This function applies a temporal difference update rule to learn the 
% successor representation (SR) and a reward estimate vector.
%
% INPUTS:
%   r      - Scalar reward at the current state
%   phi    - Feature vector of the current state (column vector)
%   n_phi  - Feature vector of the next state (column vector)
%   M      - Current Successor Representation matrix (size: [num_features x num_features])
%   R      - Current reward estimate vector (size: [num_features x 1])
%   alpha  - Learning rate for SR update
%
% OUTPUTS:
%   new_M  - Updated Successor Representation matrix
%   R      - Updated reward estimate vector
%
% PARAMETERS:
%   gamma  - Discount factor (set to 0.995, close to 1 for long-term dependencies)
%   r_alpha - Learning rate for reward update (set very small to ensure smooth updates)

gamma = 0.995;  % Discount factor for future state importance (typical: 0.995)
r_alpha = 1e-6; % Small learning rate for reward updates

%% **Update Successor Representation (SR) Matrix**
% The SR matrix M(s, s') represents the expected future state occupancy.
% The update follows a TD-style learning rule:
%
%   M(s, s') = M(s, s') + alpha * (I + γM - M) * ϕ
%
% where:
%   - `phi'` is the **current state feature vector**
%   - `n_phi'` is the **next state feature vector**
%   - `M * n_phi'` estimates future occupancy using the SR matrix
%   - `M * phi'` estimates the previous state’s occupancy
%
% This equation adjusts the SR matrix **based on the prediction error**.

new_M = M + alpha * (phi' + gamma * M * n_phi' - M * phi') * phi;

%% **Update Reward Estimate Vector** - Reward not currently used.
% The reward estimate vector `R` is updated using a temporal difference-like rule:
%
%   R(s) = R(s) + α_R * (r + γ * M * R' - M * R')
%
% This ensures `R` captures **long-term expected rewards** from the SR matrix.

R = R + r_alpha * M * phi' .* (r + gamma * dot(M * n_phi', R) - dot(M * phi', R));

end

