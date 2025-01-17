These function were originally developed for [this paper](http://www.sciencemag.org/content/337/6096/853.abstract) by JK.  However they have since been tidied up and made easier for others to use. **All these functions need to be checked thoroughly before use.**

The two main functions are [GetMaxPower2DFFTAndShuffle](#GetMaxPower2DFFTAndShuffle) and [GetPeaks2DFFT](#GetPeaks2DFFT).

----    
[GetMaxPower2DFFTAndShuffle](#GetMaxPower2DFFTAndShuffle)
| [GetPeaks2DFFT](#GetPeaks2DFFT)
| [FindPeakDirectionsInFFT2](#FindPeakDirectionsInFFT2)
| [GetFFTPowerByDirection](#GetFFTPowerByDirection)
| [GetMasked2DFFTfromRMap](#GetMasked2DFFTfromRMap) 
| [ShiftSpikePosInd](#ShiftSpikePosInd)
| [sft2_low] (#stf2_low)

----    

### <a name="GetMaxPower2DFFTAndShuffle"/> GetMaxPower2DFFTAndShuffle &#9827;  
Takes data for a single cell, and does the following...
* Build ratemap and compute the 2DFFT.   
* Find the maximum power.   
Repeat the following until the value in the final step converges:   
* Shift spike times by a random amount.  
* For the shifted spike times recompute the 2Dfft and find its peak.   
* Check to see what the 95th percentile of peak power is over all iterations so far.   

**Takes**   
Input structure has the following fields:   
`posFreq=50` frequency in Hz of position sampling   
`posXY` a matrix of size `[L x 2]` giving the (x,y) position at every moment of the trial in cm.  
`spikePosInd` a vector giving the indices into `posXY` when the cell in question fired an action potential.   
`BIN_SIZE_CM=2.5` as it says   
`MIN_SHIFT_SECONDS=20` this ensures that shuffled data is not too similar to the original   
`MAX_ITERATIONS=1000` if convergence of 95th percentile still hasn't occurred we give up   
`CONVERGENCE_THRESHOLD=0.05` maximum desired change in 95th percentile value after doing a further 100 iterations     
`FFT_SIZE=256` ratemap is zero-padded out to this size.  Powers of 2 are better, making it larger gives a smoother result.   
`MAX_FREQ=10` maximum 1d-spatial frequency (in ratemap bins), above which we ignore the power spectrum (see [GetMasked2DFFTfromRMap](#GetMasked2DFFTfromRMap)).   

**Returns**   
Output structure has the following fields:   
`maxPower` the peak average power in the unshuffled data  
`maxPowerShuffles` vector of peak average powers for all the shuffled data   


### <a name="GetPeaks2DFFT"/> GetPeaks2DFFT &#9827;  
Takes data for a single cell, and does the following...
* Build ratemap and compute the 2DFFT.   
* For each direction in the fft, find the average power.    
* Find the peak directions from these average.   
* At each of the peak directions, get a list of the 2dfft power at each wavelength.
* For each list in the above step find the peak wavelength.

**Takes**   
Input structure has the following fields:   
`posXY` a matrix of size `[L x 2]` giving the (x,y) position at every moment of the trial in cm.  
`spikePosInd` a vector giving the indices into `posXY` when the cell in question fired an action potential.   
`BIN_SIZE_CM=2.5` as it says   
`FFT_SIZE=256` ratemap is zero-padded out to this size.  Powers of 2 are better, making it larger gives a smoother result.   
`MAX_FREQ=10` maximum 1d-spatial frequency (in ratemap bins), above which we ignore the power spectrum (see [GetMasked2DFFTfromRMap](#GetMasked2DFFTfromRMap)).   
`PEAK_MAX_FRACTION=0.1`
`powerThreshold`   

**Returns**   
Output structure has the following fields:   
`pkDirections` vector of the directions with peak average power   
`half_lambda` vector of the half-wavelength of the peak spatial frequencies for the given peak directions  


### <a name="FindPeakDirectionsInFFT2"/> FindPeakDirectionsInFFT2
Takes an FFT2, computes the average power for each direction and then locates the directions with peak power.

**Takes**   
`fftPower` - the output from [GetMasked2DFFTfromRMap](#GetMasked2DFFTfromRMap), which is essentially the result of an fft2 followed by fftshift and some cropping.    
`fftDirection` -  also an output of [GetMasked2DFFTfromRMap](#GetMasked2DFFTfromRMap), which gives the direction (in degrees) of the wave vector corresponding to each bin of `fftPower`.    
`fftFreq` - same as `fftDirection`, but gives the frequency rather than direction.    
`MAX_FREQ` - used to mask a circular section of the fft (which has already been cropped to a square/rectangle).    
`thresh_coeff` - the fraction of the maximum power required for other points to be considered peaks.   

**Returns**   
`peaksDirection` direction (in degrees) of the peaks that were identified   
`peaksAvPower` the power of the identified peaks   

### <a name="GetFFTPowerByDirection"/> GetFFTPowerByDirection
Takes a list of the power in each bin of the FFT2, and a list giving the direction which each bin corresponds to, and then finds the average power for each direction.  It then circularly smooths using a guassian kernel. Note that directions should be on the interval [0 360).

**Takes**   
`fftDirection` vector of directions for the corresponding power values.  Angles should be in degrees and on the interval `[0 360]`.   
`fftPower` value from the 2dfft.   

**Returns**   
`direction` the vector `0:359`   
`avPower_smooth` vector giving the average power at the angles specified by `direction`.

### <a name="GetMasked2DFFTfromRMap"/> GetMasked2DFFTfromRMap
Performs the 2dfft on a ratemap, applies fftshift and then crops down to just the low frequencies.  Actually it uses [sft2_low](#sft2_low) rather than `fft2` plus cropping as this is faster.

**Takes**   
`rmap` the ratemap   
`S` the desired size of the fft, should be a power of 2, and bigger means more smoothing   
`maxFreq` maximum 1d-spatial frequency to calculate Fourier components for    

**Returns**   
`power_spect`  2d matrix with nans everywhere expect on the ring corresponding to the medium wavelengths.  The non-nan values give the fft2 power for the given spatial frequency.    
`direction` - gives the direction in degrees for each bin in the power spectrum   
`frequency` - gives the frequency of the wave for each bin in the power spectrum   


### <a name="ShiftSpikePosInd"/> ShiftSpikePosInd
**Takes**   
`spikePosInd` the index of pos sample corresponding to spike time
`nPos` the number of pos samples
`minShiftPos` the minimum permitted temporal shift expressed in units of pos sample period (e.g. how many 50ths of a second).

**Returns**   
`shifted_spikePosInd` the `spikePosInd` temporally shifted, wrapping the end back round to the start.



### <a name="sft2_low"/> sft2_low
If you only want the low frequencies (which we do here) it is faster to do a basic "slow" Fourier transform than to do a full two dimensional FFT.  Daniel wrote this function for that purpose and then put it on the [Matlab File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/43311-low-freq-2d-fourier-transform), where it will hopefully be checked and possibly improved by other Matlab users.
