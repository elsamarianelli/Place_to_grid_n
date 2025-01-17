function ret = power_spectrum(signal, padToPow2, binWidthSecs, freqBand,...
    maxFreq, smthKernelWidth, smthKernelSigma, s2nWdth, PLOT_ON, ymax, xmax)
% Function used by eeg_power_spectra and intrinsic_freq_autoCorr
% signal - must be mean normalised already

% #######################################################################
% [Stage 1] Get raw power spectrum

% 1a. Choose some power of 2 length to which signal shall be padded, (for fft efficiency)
nqLim = 1/binWidthSecs/2;

originalLength = length(signal);
if isnan(padToPow2)
    fftLength = 2^nextpow2(originalLength);
else
    fftLength = 2^padToPow2;
end
fftHalfLength = fftLength/2+1;

% 1b. Perform fft
fftRes = fft(signal, fftLength);

% 1c. Get power-density from fft and discard second half of spectrum
power = abs(fftRes).^2/originalLength; %dividing by length gives power density
power(fftHalfLength+1:end) = []; %for real-valued signal, power spectrum is symmetric, around central Nyquist freq, except 0 freq...
power(2:fftHalfLength-1) = power(2:fftHalfLength-1)*2; %..with half power in each side of spectrum

% 1d. Calculate freqs and crop specturm to requested range
freqs = nqLim*linspace(0,1,fftHalfLength);
freqs = freqs(freqs <= maxFreq)';
power = power(1:numel(freqs));

% #######################################################################
% [Stage 2] Smooth spectrum using gaussian kernel
binsPerHz = (fftHalfLength-1)/nqLim;
kernelLength =  round(smthKernelWidth *binsPerHz);
kernelSigma = smthKernelSigma*binsPerHz;
kern = fspecial('gaussian',[kernelLength 1],kernelSigma);
power_smoothed =  imfilter(power, kern, 'symmetric');

% #######################################################################
% [Stage 3] Calculate some metrics

%3a. Find max in theta band
spectrumMaskBand = freqs>freqBand(1) & freqs<freqBand(2);
[ret.bandMaxPower, maxBinInBand] = max(power_smoothed(spectrumMaskBand));
bandFreqs = freqs(spectrumMaskBand);
freqAtBandMaxPower = bandFreqs(maxBinInBand);

%3b. Find power in small window around peak, and divide by power
%    in rest of spectrum to get signal to noise ratio.
spectrumMaskPeak = freqs>freqAtBandMaxPower-s2nWdth/2 & freqs < freqAtBandMaxPower+s2nWdth/2;
ret.s2n = nanmean(power_smoothed(spectrumMaskPeak)) / nanmean(power_smoothed(~spectrumMaskPeak));

ret.maxFreq = freqAtBandMaxPower;
ret.power = power_smoothed;
ret.freqs = freqs;

if PLOT_ON
    PlotForDebug(power_smoothed,ret.freqs,spectrumMaskPeak,power,...
            freqBand,freqAtBandMaxPower,ret.bandMaxPower,xmax,ymax)
end


end

function PlotForDebug(power_smoothed,freqs,spectrumMaskPeak,power,...
            freqBand,freqAtBandMaxPower,bandMaxPower,xmax,ymax)
        
    if isempty(ymax)
        ymax =  min(2*max(power),max(power_smoothed));
        if ymax==0, ymax=1; end
    end
    
    % Plot "envelope" of unsmoothed spectrum, where the evelope is defiend as
    % the moving-maximum over a region 1/K times the length of the spectrum
    % (we do this because Matlab isn't good at plotting huge data series)
    K = 200;
    st = strel('rectangle',[ceil(numel(power)/K) ,1]);
    power_max = imdilate(power,st);
    isChange =   [true; power_max(1:end-1) ~= power_max(2:end)] ...
              | [power_max(1:end-1) ~= power_max(2:end) ; true];
    power_max = power_max(isChange);
    changeFreq = freqs(isChange);
    hArea = area(changeFreq,power_max);
    set(hArea,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none')
    hold on;  

    
    hArea = area(freqs(spectrumMaskPeak),power_smoothed(spectrumMaskPeak));
    set(hArea,'FaceColor',[1 0.9 0.9],'EdgeColor',[1 .5 .5])
    plot(freqs,power_smoothed,'k','LineWidth',2); 
    plot([freqBand(1) freqBand(1)],[0 ymax*10],'--b'); % times 3 is so that if you change y axis lim it'll still be ok
    plot([freqBand(2) freqBand(2)],[0 ymax*10],'--b');
    stem(freqAtBandMaxPower,bandMaxPower,'r','LineWidth',2);
    hold off;
    ylim([0 ymax]);
    xlim([0 xmax]);
    xlabel('frequency (Hz)', 'FontSize', 8);
    ylabel('power density (W/Hz)', 'FontSize', 8);
    set(gca, 'FontSize', 8)
end