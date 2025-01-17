function [direction, avPower_smooth] = GetFFTPowerByDirection(fftPower,fftDirection)
% Original code by JK.
%
% Takes a list of the power in each bin of the FFT2, and a list giving the 
% direction which each bin corresponds to, and then finds the average power
% for each direction.  It then circularly smooths using a guassian kernel.
% Note that directions should be on the interval [0 360).
%
% DM, May 2013.  Fully vectorised code and cleaned it up.
%

W = 8; % kernel width is 2W+1, i.e. 17 degrees
sigma = 13; % 13 degrees
%Note that with these param values there aren't really any tails to the filter

% prepare gaussian kernel
x = -W:W;
gausKern = 1 * exp(-x.^2/(2*sigma^2));
gausKern = gausKern./sum(gausKern);

% Bin the angle values into single-degree bins and take the mean
% of the power spectrum in each bin, i.e. average over all wavelengths
% with the same angle
directionInd = floor(fftDirection) + 1;
avPower = accumarray(directionInd,fftPower,[360 1],@mean)';
direction = 0:359; %DM, strictly speaking bin centres are probably +0.5 degrees
nDirections = numel(direction);

% smooth SF_av, padding it first so that smoothing is circular
avPower_padded = [avPower(end-W+1:end) avPower avPower(1:W)];
avPower_smooth = filter(gausKern, 1, avPower_padded);
avPower_smooth = avPower_smooth(W+(1:nDirections)); %discard the padding

end