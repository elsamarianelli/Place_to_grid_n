`LoadData` is a set of functions for loading Axona files and cut files into a Matlab structure and performing some basic post-processing to filter out spurious data in the tracking recording.  The main function is `load_trial`, see below for a full explanation, but here's an example to get you started:
```matlab 
in.flnmroot = 'experiment 25';
in.data_path = 'C:\data\rat99\';
trial = load_trial(in);
```

----   

  [load_trial](#load_trial) | [read_tetrode_file](#read_tetrode_file) | [read_eeg_file](#read_eeg_file) | [read_pos_file](#read_pos_file) | [read_binary_data](#read_binary_data) | [postprocess_pos_data](#postprocess_pos_data) | [led_speed_filter](#led_speed_filter) | [led_swap_filter](#led_swap_filter) | [read_cut_file](#read_cut_file) | [mTintHeader](#mTintHeader)

----    

### <a name="load_trial"/> load_trial &#9827;     
See above for an example of basic usage.

**Takes**      
Input structure has the following fields:    
`data_path` e.g. 'C:\example\', i.e. with trailing slash    
`flnmroot` e.g. `'R1884 2008-01-27 3'`, i.e. whatever you saved your set file as minus the file extension    
`GET_TRODES=[1 2 3 4]` list of which tetrodes to load data for, can be an empty array. See `spikeData` below.    
`cutSuffix='_'` this is the bit that comes between the experiment name and the `x.cut`, where `x` is the tetrode number.    
`USE_CLU_FILES=false` if true cut file names end in "clu.#tet" rather than "#tet.cut". Note that you will probably need to set the `cutSuffix` to `''` rather than the default of `'_'`.    
`GET_POS=true` if false pos file data is not loaded    
`GET_WAVEFORMS=false` if true waveforms as well as timestamps are loaded form tetrode files    
`GET_CUT=true` if false cut file is not read    
`GET_EEG=false` if true all eeg data is read (except high frequency eeg files)   
`POS_PROCESS_MAX_SPEED=5` meters per second, see [postprocess_pos_data](#postprocess_pos_data).    
`POS_PROCESS_BOX_CAR=0.2` again see [postprocess_pos_data](#postprocess_pos_data).    
`POS_HEADER_MODIFY` a list of key-name, key-value pairs with which to modify the header for the pos file at the point it is loaded.    
`ACC_BOX_CAR=5` the width of the acceleration box car smoothing in bins.    
`CONVERT_EEG_TO_VOLTS=true` multiplies eeg values by a scale factor see [read_eeg_file](#read_eeg_file) for details.    
`CONVERT_POS_TO_CM=true` when set to true, divides pos xy values by `pixels_per_meter` and multiplies by `100`.  Note that the unit of length for speed and acceleration is always centimeters not pixels.   
`DISCARD_MINUTES` an `[n x 2]` matrix giving a list of `n` start and end times in minutes for which to discard data.  If you give a `NaN` in the second column only the minute in the start column is discarded, i.e. for input `[2 5; 8 NaN; 12 13]` minutes `2`, `3`, `4`, `8`, `12`, and `13` will be discarded.  All pos, tetrode and eeg data is discarded for the given minutes.     
`GET_EEG_RANGE=false` when set to true, the EEG mean and std are calculated and combined to give the range: `range=mean + 2*std`, otherwise the range is not set.     


**Returns**    
Output structure has the following fields:    
`setfile_header` an mTintHeader object - see [mTintHeader](#mTintHeader).    
`eegData` a structure containing headers and data for each standard-rate eeg channel    
`trialData` a structure giving the trial start and end time in seconds    
`posData` a structure containing the pos file header and the raw LED data as well as the post-processed data which includes position, direction, speed, and acceleration.    
`spikeData` and _array_ of structures, one for each tetrode (the array is of size 4 unless there is a larger value in `inStruct.GET_TRODES`). Note that even if you only ask for data from tetrode 3, say, you will still get an array of size 4, with the third element having the data in question. It is designed this way to ensure the user is aware which tetrode they are working on. The structure contains the tetrode header, spike timestamps and waveforms as well as the corresponding index in the pos data. The cut data for the tetrode is also given in the structure.    
   
### <a name="read_tetrode_file"/> read_tetrode_file    
Uses [read_binary_data](#read_binary_data) to read a tetrode file and then converts the array of bytes into a useful Matlab structure.    

Two forms:    
`[header, timestamp, waveforms] = read_tetrode_file(flnm)` does read waveforms    
`[header, timestamp] = read_tetrode_file(flnm)` does not read waveforms    
   
**Takes**    
`flnm` - fully qualified filename including path plus extension e.g. '.1', '.2'    

**Returns**    
`header` - header information in an [mTintHeader](#mTintHeader) object.    
`timestamp` - timings for each spike. Array of size `[num_spikes x num_chans]`, though all columns are identical (since they are obtained with `repmat` applied to one column of data).    
`spike` - spike waveforms for each channel `[num_spikes x samples_per_spike x num_chans]`. If you do not want the spike waveforms do not specify a third output, this will significantly improve the speed of the function.    

### <a name="read_eeg_file"/> read_eeg_file    
Uses [read_binary_data](#read_binary_data) to read an eeg file and then uses the value of sample_rate in the header to generate timestamps. Usually each sample is a single byte, but it can be two. If requested it will convert the raw values to voltages.    
`[eeg, header, scaleFactor]=read_eeg_file(flnm[,eegNum,setfile_header,CONVERT_TO_VOLTS]))`    

**Takes**    
`flnm` \- the full filename    
The other three inputs are only required if you want to convert the raw eeg to volts.    
`eegNum` \- which of the eeg channels this file corresponds to (normally 1,2,3, or 4)    
`setfile_header` \- the header from the set file    
`CONVERT_TO_VOLTS` \- conversion to volts is performed if this flag is true.    

**Returns**    
`eeg` \- array of size `[num_samps x 2]` where `eeg(:,1)` gives the eeg values and `eeg(:,2)` gives the timestamps.    
`header` \- header information in an [mTintHeader](#mTintHeader) object. Normally contains (among other things) number of channels and sample rate.    
`scaleFactor` \- equal to `1`, unless a conversion to volts was performed, in which case it gives the factor used.    

**Notes**    
The eeg sample rate is normally 5 times the tracking rate, i.e. only every fifth eeg sample gets a tracking data packet.    
The conversion to volts is done using a factor composed of three separate numbers:    
* The `'ADC_fullscale_mv'` key in the setfile header divided by `1000`.    
* The `'gain_ch_X'` key in the setfile header, where `X` is the actual channel number minus one, not the eeg channel number. To get the actual channel number from `eegNum` the `'EEG_ch_Y'` key in the setfile header is used, where `Y` is `eegNum`.    
* The maximum int value, this depends on the `'bytes_per_sample'` key in the eeg header. Values are normally stored as signed ints of either 1 or 2 bytes, and must be scaled to the interval [-1,1] to be used with the other two factors.    

###  <a name="read_pos_file"/>  read_pos_file    
Uses [read_binary_data](#read_binary_data) to read a pos file and then converts the array of bytes into a useful Matlab structure. Does not perform any post-processing, except to remove unused LEDs (based on having only zeros or NaNs in the relevant column).    
`[led_pos, led_pix, header] = read_pos_file(flnm)`    

**Takes**    
`flnm` \- the full file name of the pos file    

**Returns**    
`led_pos` \- an array giving the locations of the LEDs at each tracking event. For `n` tracking events, `L` LEDs, it will be of size `[nxLx2]`, where the `2` is for the two spatial dimensions.    
`led_pix`    
`header` \- header information in an [mTintHeader](#mTintHeader) object. Normally contains (among other things) pixels per meter and window max/min.    

**Note**    
Unlike with eeg loading, conversion to meaningful units (in this case centimeters) is not an option within this function, instead it is dealt with in [load_trial](#load_trial).    

### <a name="read_binary_data"/> read_binary_data    
Opens a file, finds the section before the letters "data_start" and parses it into a header variable. The contents of the file between the letters "data_start" and "data_end" is returned as a big array of single bytes.    
`[header, RawBinaryData] = read_binary_data(filename,type)`    

**Takes**    
`filename` \- the full filename to read in    
`type` \- one of: `"set", "pos", "eeg", "tet", "tet_times"`    

**Returns**    
`header` \- header information in an [mTintHeader](#mTintHeader) object.    
`RawBinaryData` \- an array of bytes    

**Notes**    
The function searches for "data_start" from the beginning of the file and "data_end" from the end of the file, going backwards (well, in batches going forwards starting from the end).    
The header format in the file is [key_name] [space] [key_value] [newline].    
If `type==".set"` the file only contains a header (not even a "data_start" marker) so the whole file is parsed in the header variable.    
DM. Previously the function read all of the file into memory and then passed it to a helper function to find the data start and data end markers, whereas now it uses `fseek` to read straight into the target memory (you may recognize some of the code from a recent version of [read_tetrode_file](#read_tetrode_file) - indeed [read_tetrode_file](#read_tetrode_file) now uses this function for loading data rather than doing it itself)    

### <a name="postprocess_pos_data"/> postprocess_pos_data    
Performs the following steps to calculate position at each timestep (normally 0.02s):    
* call and then apply [led_swap_filter](#led_swap_filter)    
* apply [led_speed_filter](#led_speed_filter)    
* interpolate values in filtered out regions of data    
* apply box-car smoothing (moving average)    
Finally, use the position values to calculate speed and direction.    
`[xy, dir, speed, times, jumpyPercent, n_leds, dir_disp] = postprocess_pos_data(posdata, max_speed, box_car, setfile_header)`    

   
**Takes**    
`posdata` \- structure with the following fields:    
&nbsp;&nbsp;&nbsp;`.led_pos` \- an array giving the locations of the LEDs at each tracking event. For `n` tracking events, `L` LEDs, it will be of size `[nxLx2]`, where the `2` is for the two spatial dimensions.    
&nbsp;&nbsp;&nbsp;`.led_pix` \- which (of the multiple) LEDs to apply the filtering to. Number corresponds to the second dimension in `led_pos`    
&nbsp;&nbsp;&nbsp;`.header` \- pos file header    
`max_speed` \- max speed in meters per second used for [led_speed_filter](#led_speed_filter).    
`box_car` \- width in seconds to use as a box-car smoothing kernel (i.e. width of moving average). Give 0 for no smoothing.    
`setfile_header` \- the main header for the trial, given in the .set file. (Note the additional header for the tracking data passed in `posdata.header` field above.)    

**Returns**    
`xy` \- position (in pixels - ij format),    
`dir` \- direction in degrees anticlock wisefrom x axis,    
`speed` \- in cm/s (using pixels per metre and position sample rate in `posdata.header`).    
`times` \- in s, vector of time of each pos sample starts with e.g. `[1/50, 2/50, ...]`    
`jumpyPercent` \- percent of pos points that were excluded because of too high speed    
`n_leds` \- number of LEDs used during recording    
`dir_disp` \- direction (degs anti-clock from x-axis) derived from heading. In case of `n_leds == 1` then `dir_disp==di`    

   
### <a name="led_speed_filter"/> led_speed_filter    

Starting with tracking samples `a,b,c,d,...` the function calculates the running speed required to get from `a` to `b` given the time between the two observations. Where this speed is greater than some chosen threshold, it is assumed that the sudden jump is a tracking error rather than a true sprint by the animal. Tracking sample `b` is set to `NaN` and the algorithm goes on by testing the speed required to get from `a` to `c`. Thus the resulting list of tracking samples should have any spurious jumps filtered out and set to `NaN`.    

`[n_jumpy, led_pos] = led_speed_filter(led_pos, max_pix_per_sample, led)`    

**Takes**    
`max_pix_per_sample` \- speed threshold, defined as maximum pixels moved between tracking events (usually 0.02 seconds).    
`led_pos` \- an array giving the locations of the LEDs at each tracking event. For `n` tracking events, `L` LEDs, it will be of size `[nxLx2]`, where the `2` is for the two spatial dimensions.    
`led` \- which (of the multiple) LEDs to apply the filtering to. Number corresponds to the second dimension in `led_pos`.    

**Returns**    
`n_jumpy` \- the number of tracking events set to `NaN`    
`led_pos` \- the filtered version of the input `led_pos` array    

**Notes**    
DM - have added a complicated optimized version of the function which is heavily commented but not exactly comprehensible. The slow version remains in the file, but is commented out.    

### <a name="led_swap_filter"/> led_swap_filter    

Checks for instances of two LEDs swapping or big one replacing little one when the big one gets obscured.    
The algorithm looks at the distribution of the number of pixels tracked for each LED over the trial as well as calculating the distance between successive points.    

`swap_list = led_swap_filter(led_pos, led_pix)`    

**Takes**    
`led_pos` \- an array giving the locations of the LEDs at each tracking event. For `n` tracking events, `L` LEDs, it will be of size`[nxLx2]`, where the `2` is for the two spatial dimensions.    
`led_pix` \- an array giving the number of pixels tracked for each LEDs at each tracking event. For `n` tracking events, `L` LEDs, it will be of size`[nxL]`.    

**Returns**    
`swap_list` \- list of tracking events where swap is thought to have occurred.    

   
### <a name="read_cut_file"/> read_cut_file    

`[exact, unique_clusters, n_user_clusters] = read_cut_file(filename)`    

**Takes**    
`filename` \- full file name of .cut file    

**Returns**    
`exact` \- vector of cluster group ids, one number for each event on the tetrode. Cluster group ids are in the range 0\-`n_user_clusters`.    
`unique_clusters` \- list of non-empty clusters (including 0).    
`n_user_clusters` \- the largest cluster group id (less than 30, unless you've written your own cutting tool\!).    

### <a name="mTintHeader"/>  mTintHeader    
This is a class not a function.    

**Constructor**    
To create an `mTintHeader` object from a file, find the header text and then do:    
`myHeader = mTintHeader(header_text_from_file,[modifications])`    
`modifications` is an optional cell array of key-name key-value pairs which will be used to modify keys (Matching is exact, if the key name is not found a new key is added.)    

**Methods**    
You can then use `myHeader.KeyValue(key_name,[format])` to get a key value.   
To match a key name exactly, rather than just matching on the start of the key name, use `myHeader.KeyValueExact(key_name,[format])`.    
`[format]` is optional, by default a cell array of string is returned. Alternatively:    
&nbsp;&nbsp;`'num'` \- convert key value string(s) to a double    
&nbsp;&nbsp;`'string'` \- return key(s) as a single space-padded char array rather than a cell array    
&nbsp;&nbsp;`'int8', 'uint32', 'single', ...` \- cast numeric output to the given Matlab type    
&nbsp;&nbsp;`'echo'` \- displays the values on command line.    

To view the whole header do:    
&nbsp;&nbsp;`myHeader.ToString`    

You can combine `ToString` with the `modifications` input in the constructor to make a new header with modified keys:    
&nbsp;&nbsp;`newHeader = mTintHeader(oldHeader.ToString,modifications)`    
You could also potentially write out the modified header string to file. However it is advisable to write to a new copy of the file rather than overwritting the original.    

To test if one or more keys are present do:    
&nbsp;&nbsp;`isPreesnt = myHeader.KeysPresent({'one_thing','another_thing'})`    
The output is a vector of logicals indicating whether an exact match was found for each of the key names.    

To get multiple key values in one go do:    
&nbsp;&nbsp;`[one_value,another_value,...] = myHeader.KeyValueMultiple('one_thing','another_thing',...,'num')`    
Matching is exact and you must specify the format using the final argument, i.e. it is not optional here (it can't be optional because the function wouldn't be able to tell what you meant).
