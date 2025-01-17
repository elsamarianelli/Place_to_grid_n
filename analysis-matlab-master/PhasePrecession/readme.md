### Phase Precession Documentation

1.	[Download](https://github.com/UCL/mTint/archive/master.zip) code package from the UCL/mTintGitHub repository.
2.	Unpack and add all folders/sub folders to your MATLAB path.
3.	Separately, you will need the data you want to analyse and a tab delimited text file that contains columns of relevant metadata, with each row corresponding to a single cell. Examples are provided.
4.	Before you can get started you will need to create two functions that are relevant to your data. One is to create an input array of structures that has one element per cell. As an example, see `createInputArrayForSargoliniDataSet` or `createInputArrayForOriginalGridCellDataSet`. At a minimum, your input array must contain a field `type`, with entry `'g'` or `'p'`, to indicate grid cells or place cells, along with all the information you need in the second function.
5.	The second function is the loading function relevant to loading your data. (See the functions `loadSargolini` and `loadOriginal` as examples). Your load function takes two inputs, one is the array created in (4), the other is an FIR filter created in `ppSettings` (see `SETTINGS.LOADOPTIONS` in `ppSettings`) to filter the EEG signal. This loading function is required to output the following:
  *	`xy` - the position data in xy co-ordinates (preferably in cm)
  *	`dir` - the direction data preferably calculated from the position (not 2-spot direction)
  *	`spd` - the speed preferably calculated from position as well.
  *	`eeg` - the raw EEG signal
  *	`phase` - the EEG phase (requires the smoothed EEG below)
  *	`filteredEEG` - the smoothed EEG signal (this is where the FIR filter comes in)
  *	`psr` - the position sampling rate (for axona, obtained from the pos header)
  *	`esr` - the EEG sampling rate (for axona, obtained from the eeg header)
  *	`spikeTS` - the timestamps of spikes fired by a cell
  *	`type` - the type of cell eg, `'g'` or `'p'`.
6.	Once these are ready, you need to go through `ppSettings` and make sure they are correct. For the most part, you can leave them as they are. If you want to play with the settings, you'll need to study the relevant functions carefully to understand what's going on. But many things are fairly straightforward and there is some commenting to help you along and reasonably sensibly named variables. In `ppSettings`, you have to set `SETTINGS.outPath`. This is where all the output gets saved. The way the system works is that it saves the data one cell at a time to avoid running out of memory, then at the end combines all the files into one.
7.	You are now ready to run the main function `ppMain`, see comments within the file to see how to use.
