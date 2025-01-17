function [power_spect, varargout] = GetMasked2DFFTfromRMap(rmap, S, maxFreq)
%
% Original code by JK.
%
% DM, May 2013.  Cleaned up.
%

% mean normalise and remove nans
nanmean_rmp = mean(rmap(~isnan(rmap)));
rmap = rmap-nanmean_rmp;
rmap(isnan(rmap)) = 0;

% compute the 2d-FT for frequencies 0 to maxFreq
if nargout == 1
    power_spect = sft2_low(rmap,maxFreq,S);
else
    [power_spect, freqs_1,freqs_2] = sft2_low(rmap,maxFreq,S,1);   
    freqs_1 = fftshift(freqs_1);
    freqs_2 = fftshift(freqs_2);
    varargout{1} = mod( bsxfun(@atan2d, freqs_1, freqs_2'), 360); % direction in degrees
    varargout{2} = sqrt(bsxfun(@plus, freqs_1.^2, freqs_2'.^2)); % 1d freq
    
end
power_spect = fftshift(power_spect.*conj(power_spect)); 


end
