function f= fftVars
%FFTVARS Structure of common vars used by family of functions doing 2d FFT


% --- BASIC --------------------------------------------------------------
f.binSizeCm=2; %Size of spatial bins in cm - converted to pix using ppm
f.minShuf=20; %Min amount to shuffle spikes by in sec [JK used 20]


% --- CREATE FFT ---------------------------------------------------------
%FFT is zero paded to be square with sides of this length. Ideally for grid
%cells need something in the order of 2048 to 4096 to ensure accuity of at
%least 1cm in 60cm range for 2cm bins. JK used 256.
% Note changing f.fftSize requires changes to be made to the following:
% f.smth2dSig [if used], f.maxRad, f.thBin [in direct effect], f.rSig
% [smoothing sigma in FFT freq components]
f.fftSize=512; %[1024] 1024 is sufficient and probably 512, 256 is too low

%Smth 2D FFT power spc with guassian of this sig (0=no smth). Note best not
%to smooth in 2D coords
f.smth2dSig=0; 

%Can return only the central part of the power spec (after low freq shifted
%to centre). If so maxRad defines how many bins from the centre are
%returned. JK used 40bins but that was for 2.5cm spatial bins and 256N FFT
%(this excluded wavelengths less than ~16cm). For 2cm spatial bins and
%2048N to exclude wavelengths less than 10cm need 410bins radius.
%NB makes big difference to result if radius excludes some of the main
%components and of course the postion of these components moves when
%fftSize is changed
%For 256N 16cm ~40
%For 1024N 16cm ~140
%For 2048 16cm ~200
%NB have changed now define minWaveCm i.e. the minium wavelength in cm that
%should be returned, this is converted to a maxRad based on the binSize and
%size of FFT.
f.minWaveCm=12; %JK value was equal to about 14cm
f.maxRad=ceil((1/f.minWaveCm) * (f.binSizeCm * f.fftSize)); %[410]

% --- POLAR FFT ----------------------------------------------------------
%Convert 2D FFT to polar and then smooth in polar space (polar FFT is an mn
%mat where dim m is angle in degs measured anticlock from x-axis and n is
%FFT component so similar for 2D FFT. Note polar FFT is truncated at
%f.maxRad defined above.
f.thBin=1; %Ang bin in deg [1] JK used 1.  <1 requires fftSize>256
f.rBin=1; %Radial bin size - best to leave as 1 to match 2dFFT
%NL. Smoothing. JK only smoothed in angluar space with very wide guassian
%(sigma = 13deg) so effectivly a box car. I might get away with smaller as
%doing in 2D. Note my value is in degs so takes account of f.thBin
f.thSig=10; %Sigma of gaussian for ang smoothing [8] .NB in degs
% Radial sigma - JK didn't use but my value is in FFT components so adjusts
% for f.rBin but will need to be changed for larger f.fftSize
f.rSig=4; %Radial sigma in FFT components [3]


% --- SHUFFLING & THRESHOLD ----------------------------------------------
%Number of iterations used to get shuffled dist [JK used 1000 but with
%stopping criterion]
f.nIt=1000; %1000 is best 500 OK
%Percentile threshold to use for significant FFT components in 2d FFT. 
%This percentile value is selecte from shuffled distribution max power. 
%Note JK uses 50 but I was getting better results with 95.
f.pThresh=50;
%Secondary peak threshold. Once collapsed to angular FFT need to find peaks
%but only return peaks which are above some proportion of the max peaks. JK
%used 0.1 (i.e. 10%)
f.sndPeakCoef=0.2; %[0.1]
%JK also required the peaks to be a certain number of bins distant from
%each other (she used 10deg). Note my value is in degs so takes acount of
%f.thBin
f.minPeakDist=10; %Min seperation of peaks in deg

end

