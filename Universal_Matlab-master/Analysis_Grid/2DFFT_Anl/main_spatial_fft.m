function [ ps2D, psPol, fCompAngDeg,fCompWaveCm ] = main_spatial_fft( posData, spkInd )
%MAIN_SPATIAL_FFT Bins, gets 2D FFT, surrogates & sig components
% Master function for analysing FFT components of spatial cells - calls
% other functions to do most of the work. Starts with raw positional data
% and spike indicies, then bins and gets 2D FFT. Also shuffles spikes and
% repeates process to get confidence intervals.
%
% ARGS
% posData - pos branch of tint global data structure produces by mtint read
%           in functions i.e. tint.pos
% spkInd - [nSpkx1] index that spikes occured at where values are index
%           into posXY. e.g. data.tetrode.pos_sample
% Also requires access to function fftVars which carries a variable struct
%
% RETURNS
% ps2D - [N x N] 2D FFT of ratemap thresholded at 50% conf interval (low
%       freq at centre)
% psPol - [nDeg x maxRad N] ps2D transformed into radial coords (low freq
%       on left)
% fCompAngDeg - [nComponents x 1] angle in deg (antiClock from x-axis) of
%       sig FFT components. Note these lie at 90deg to grid axis.
% fCompWaveCm - [nCompoents x 1] wavelength in cm of FFT components will be
%       shorter than grid scale (for reg grid
%       fCompWaveCm=gridScale*cosd(30)


% ------------------------------------------------------------------------
% --- HOUSE KEEPING ------------------------------------------------------
% ------------------------------------------------------------------------
goPar; %Start parallel instances - used to speed shuffling
f=fftVars; %Load vars structure from this folder

nPosPt=length(posData.xy);
ppm=str2double(posData.header{strcmp('pixels_per_metre', posData.header(:,1)),2});
spBin=ppm/100*f.binSizeCm; %Spatial bin size in pix for a fixed size in cm
clear ppm

% ------------------------------------------------------------------------
% --- MAIN ---------------------------------------------------------------
% ------------------------------------------------------------------------

%1)Start by binning data to get a spatial ratemap (unsmoothed).
binXY=single(bin_pos_data('position', spBin, posData, 1:nPosPt));
binSpk=single(bin_pos_data('position', spBin, posData, spkInd));
rm=binSpk./binXY; %ratemap with no smoothing


%2)Now get the real power spectra
[ ps2D ] = get_2d_fft( rm, f.fftSize, f.smth2dSig, f.maxRad );


%3)Make surrogates to compare against real power spec.
shufMaxPower=zeros(f.nIt,1); %Preallocate- max from each shuff
% tmpAll=zeros(f.fftSize, f.fftSize, f.nIt); %Store all ps  - for debug
minShuf=f.minShuf*50; %Min shift in pos pt from seconds

parfor nn=1:f.nIt
    %First shift spkind by some random amount, make new rm & power spec
    shufSpkInd=sf_shuf_spk_ind(spkInd, nPosPt, minShuf);
    shufBinSpk=single(bin_pos_data('position', spBin, posData, shufSpkInd));
    shufRm=shufBinSpk./binXY; %shuffled ratemap with no smoothing
    shufPs2D  = get_2d_fft( shufRm, f.fftSize, f.smth2dSig, f.maxRad ); %shuf ps
    
    shufMaxPower(nn)=max(shufPs2D(:)); %Store max power from each ps
%     tmpAll(:,:,nn)=shufPs2D; %Store all power spec - for debug    
    
end
%Get two thresholds 50th and 95th
crit95=prctile(shufMaxPower,95);
% crit50=prctile(shufMaxPower,50);%Get 50th percentile
clear shuf* minShuf

%4) Apply critical values to real power spec and remove zeros
tmp=ps2D-crit95;
if ~any(tmp>0) %Check if there are peaks above 95th perc if so continue
 [ps2D, psPol, fCompAngDeg, fCompWaveCm]=deal([]); %Set to all empty and return
 return
end
ps2D=ps2D-crit95; %Threhold at 50th percentile
ps2D(ps2D<0)=0;
clear tmp crit95 crit50

%Now get polar fft and smooth it - set everything to single for speed
tmpThSig=f.thSig/f.thBin; %Smoothing in angular space in bins from degs
tmpRSig=f.rSig/f.rBin; %Smoothing in freq space in bins from FFT components
[ psPol, th, frq ] = ...
    get_polar_fft( ps2D, f.maxRad , f.thBin, f.rBin, [tmpRSig, tmpThSig] );
clear tmp*

%Collapse polar rm to just mean power vs angular and get the peaks in deg
angMean=mean(psPol,2); %Colapse to angluar space - smoothed in previous step
sndPeakThresh=max(angMean).*f.sndPeakCoef;
[angPeakPow,pInd]=... %angPeakPow is not retured - value at peak in ang plot
    findpeaks(double(angMean), 'minpeakheight', sndPeakThresh, 'minpeakdistance',f.minPeakDist/f.thBin);
fCompAngDeg=th(pInd); %Angle of fourier comp in deg

%Then from polar FFT find for each angle identified the freq componenet
%with the highest power hence get the freq and wavelength. Note wavelength
%will be less than for the grid scale for a totally regular grids
%wavelenght of fourier component will be cosd(30)* grid scale i.e. 0.866
[polPeakVal, peakFComp]=max(psPol(pInd,:),[],2); %polPeakVal not returned - val at peak in 2d
peakFComp=peakFComp-1; %Subtract 1 becuase first is 0th component
fCompWaveCm=f.fftSize./((1/f.binSizeCm).*peakFComp); %Get wavelgh of each FFT comp

%But note grid axis will 


%Image for debug
% subplot(1,3,1), imagesc(ps2D), axis square;
% subplot(1,3,2), imagesc(psPol);
% subplot(1,3,3), plot(th,angMean);





end

function spkInd=sf_shuf_spk_ind(spkInd, nPosPt, minShuf)
%Shift all spk ind by some random value not less that +/- some value then
%wrap to be meaningful values
% ARGS
% spkInd - index of spikes in pos points
% nPosPt - total number of pos points - max ind cannot exceed htis
% minShuf - min shuffle allowed in pos pts +/-

randShft=ceil(rand(1).*(nPosPt-minShuf)) + minShuf; % between min shuff and maxind
spkInd=mod(spkInd+randShft, nPosPt); %Add rand shift and wrap
spkInd=max(1,spkInd); %Remve zeros and set to1.
end

