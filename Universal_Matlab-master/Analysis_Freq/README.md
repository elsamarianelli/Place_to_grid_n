#Analysis_Freq

Code for the analysis of oscillatory signals (e.g. LFP, theta-band modulation of spiking etc.). 

Note for the units of LFP to meaningful it needs to be converted into volts but for single to noise (s2n) measures that absolute units don’t actually matter. Also note that Axona recording systems default to saving standard EEG files at 8bit 250hz (.eeg) but can also save a higher precision higher rate file (.egf at 16bit 4.8kHz), these require different functions to read in.


## List of contents
Might not be complete

**DETECT_RIPPLES** Very basic script for detecting the presence of ripples in the LFP. Basically takes the LFP (because ripples are high frequency ~180Hz requires the higher frequency LFP files i.e. *.egf). Filters in the ripple band, calculates the analytic signal (Hilbert), finds instantaneous power and then then looks for periods that are some number of standard deviations above mean. Note ripples can be confused with chewing artefacts. Calls *eeg_instant_freq_power.m* to calculate that analytic signal.

**EEG_CONVERT_TO_VOLTS** EEG Axona files are stored as either an 8bit signed integer (250Hz standard eeg) or 16bit integer (4.8kHz high frequency eeg files). In each case these need to be converted to volts if you want to use meaningful units for your LFP. The exact conversion factor has changed between different Axona systems but has remained constant for all the USB based systems. This code checks the header information (saved in the .set file) to determine what conversion factor to use for the conversion. The other factor that has a bearing on the conversion to volts is the gain applied by the user to the LFP channel - again this is specified in the .set file. Hence the function requires the eeg time series, the channel the eeg was recorded on (so it can look up the gain), the set file header, and if this is a 16bit eeg file a fourth argument ‘16bit’ also needs to be passed in otherwise the function assumes this is 8bit.

**EEG_INSTANT_FREQ_POWER** Calculates the instantaneous frequency and power for an LFP signal using the Hilbert transform to get the analytic signal. User specified parameters allow the signal to be filtered around a specific band. Note Hilbert can be very sensitive to filtering in the correct band - if signal actually lies outside of the band or even towards one edge then results can be biased.

**EEG_POWER_SPECTRA** Basic function for analysis of the LFP and construction of a power spectra from the same. Analysis can be limited to specific sections of the LFP (by specifying which pos samples to use), also returns that signal to noise (s2n) which is an expression of how much power lies in a particular band vs the rest of the band analysed (e.g. if s2n is returned for the theta band [6,12]hz then the peak in that band is found, a region round is detected and the power contained therein normalised by the power in the rest of the band [2,125]Hz. NB the upper limit of the frequency band analysed is specified by the user and this will affect the s2n since there isn’t a great deal of power in higher bands).

**GETEGF** Used to read in the high frequency and higher resolution *.egf files saved by DACQ i.e. the standard LFP files are saved at 250Hz and 8bit, the *.egf files are saved at 4.8kHz and 16bit. In general the higher resolution is always useful but is critical for analysis of high frequency signals (such as ripples).



