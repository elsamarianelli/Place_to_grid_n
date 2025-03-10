function [ power, freqs, relFreq, relPower ] = circ_FTT( dirRate, smthKern )
%CIRC_FTT Identifies freq present in polar data i.e. HDC tuning
% Code primarily designed to detect multimodal tuning in head direction
% tuning but might be useful for other analyses too. Take HDC data - which
% should be unsmoothed - and gets FFT. FFT naturally deals with fact that
% data is circular since this is an asumption of the FFT anyway.
% 
%
% TAKES
% dirRate       directional firing rate should be unsmoothed
%
% smthKern      [1] description of kernel to smooth data with. Leave empty
%               for no smoothing. Otherwise specified width of gaussian in
%               Hz that is applied to the power spectra. Either leave empty
%               [] or try 0.25.
%
% RETURNS
% power         Power represented at each frequncy component - arbitary
%               units. Note signal is mean normalised so 0 freq should have
%               0 power 
%
% freqs         Frequencies coresponding to the values in power. To display
%               data do 'plot(freqs,power)'
%
% relFreq       Frequncies that are compared against 6Hz
%
% relPower      For each of the frequncies specified in relFreq (these are
%               hard coded on line 71) the power in those components
%               relative to the power in the 6Hz components (i.e. 6Hz
%               component is set to 1 and the others are expressed relative
%               to that - we hope to see relative power <1 in these other
%               components)



% --- HOUSE KEEPING ------------------------------------------------------
dirRate     =dirRate(:); %col vector
sampFreq    =length(dirRate); %number of bins is sampling freq


% --- MAIN ---------------------------------------------------------------

% ----
%First do the FFT and get matching frequency components
%Mean normalise - we don't need the DC component
dirRate     =dirRate - mean(dirRate);

%Actually do the FFT
% fftLength   =2^nextpow2(sampFreq); %FFT is faster when n is power of 2
%dirRate is padded to 360 bins (if it isn't already) - hence meaning a freq
%of exactly 6Hz can be detected
fftLength   =360;
fftData     =fft(dirRate, fftLength); %result is complex

%Calculate power per component - we don't need phase
power       =(fftData.*conj(fftData))/length(fftData);% eqiv of taking the square modulus
power       =power(1: fftLength/2 +1); %Start (zero comp) to middle n/2 +1

%And the matching frequency components
freqs       =sampFreq/2  *  linspace(0,1,(fftLength/2)+1);


% ----
%Second smooth power spectrum if required
if ~isempty(smthKern)
    kernSig     =round(smthKern * ((fftLength/2 +1)/(sampFreq/2))); %Sigma in bins - convert from Hz
    kern        =fspecial('gaussian',[kernSig*5,1], kernSig);
    power       =imfilter(power, kern, 'symmetric');
end


% ----
%Third calculate relative power at other frequencies that are candidates
%for modulo directional modulation (i.e. 2, 3, 4, 5, 7
relFreq         =[2,3,4,5,7]; %Frequencies to compare against
relPower        =power./power(7); %Set power at 6Hz to 1
relPower        =relPower(relFreq+1);


end

