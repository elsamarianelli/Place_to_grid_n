function [bits_per_spike,bits_per_sec] = skaggsInfo(rates, times)
% Returns Skaggs et al's estimate of spatial information.
%
% Code to return the spatial information conveyed by a neuron. Typical
% measure is bits_per_spike but bits_per_sec is sometimes used. Note should
% really be used with adaptive binning to ensure a smooth ratemap otherwise
% values can be spurious.
%
%   I = sum_x p(x) r(x) log(r(x)/r)
%
% then divide by mean rate over bins to get bits per spike.
%
% Binning could be over any single spatial variable (e.g. location, direction, speed).
% NB this version adapted from Ali and Tom. Tested and agreed on.

if max(max(rates)) <= 0
    bits_per_spike = NaN;
    bits_per_sec = NaN;
    return
end

if (size( rates ) ~= size(times))
    error('Size mismatch between spikes and position inputs to skaggs_info');
end

rates(isnan(rates)) = 0;   % remove background
times(isnan(times)) = 0;

if(size( rates ) > 1)
    % turn arrays into column vectors
    rates = reshape( rates, numel(rates), 1);
    times = reshape( times, numel(times), 1);
end

duration = sum(times);
mean_rate = sum(rates.*times)./duration;

p_x = times./duration;
p_r = rates./mean_rate;
dum = p_x.*rates;
ind = find( dum > 0 );
bits_per_sec = sum(dum(ind).*log2(p_r(ind)));   % sum( p_pos .* rates .* log2(p_rates) )
bits_per_spike = bits_per_sec/mean_rate;