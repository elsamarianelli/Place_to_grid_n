### fix_binary_data
If DACQ crashes during a recording and it's not possible to load the temporary files in Tint/MATLAB it may be possible to fix the files.  Fixing generally means adding on a "data_end" token to the end of the data and putting in nSpikes/nPos values into the header.

Other fixes may be possible/neccessary, but currently this just works on 'tet', 'pos' and 'set' files.

In addition to writing out a corrected file, the function returns values in the same way as [read_binary_data](../LoadData#read_binary_data).

Example:   
```matlab
[header, RawBinaryData] = fix_binary_data(['C:\my data\failed trial.3','tet','.bad',1200)
```

**Takes**
```matlab
[header, RawBinaryData] = fix_binary_data(filename,type,suffixForOld,duration)
```

`filename` - full path name of a single DACQ file, i.e. tet or pos file (including extension)     
`type` - currently supports `'tet'` `'pos'` and `'set'` only.   
`suffixForOld` - if file is corrected, the origial will be copied and saved with this string append to the name e.g. `'.bad'`. (Note that if such a file already exists it will not be re-coppied.)   
`duration` - if provided it will set this value in the header and (in the case of pos) extend/shirnk the data to match.   

### merge_binary_data
If you had to split one trial into two or more parts, but wish to treat them as a single trial, you can use this function to join the subtrials together and output a new set of set/pos/tet/eeg files consisting of all the data.
Most of the header values are taken from the first of the subtrials, with a few (e.g. `duration`, `num_spikes`) set to the total value from all the subtrials.  Some additional header values are added, with `_original` suffix, listing all the original values in the subtrials, separated by askterisks.  And a special key `merged_nFiles` is added which states how may subtrials were combined to make this file.  Note that if you combine trials in multiple stages this will only show the number joined at the last stage (maybe it ought to give `total + nSubtrials` from the subtrials..but it doesn't).

Example:
```matlab
flnm = 'C:\somehwere\2141 trial3';
merge_binary_data([flnm '.pos'],{[flnm 'a.pos'],[flnm 'b.pos']},'pos','.old');
```

**Takes**
```matlab
[header, RawBinaryData] = merge_binary_data(destFilename,srcFilenames,type,suffixForOld)
```

`destFilename` - full path, name and extension for output file.   
`srcFilenames` - cell array of strings giving full path name and extensions of subtrials to be merged.    
`type` - currently supports `'tet'` `'set'` `'eeg'` and `'pos'`.     
`suffixForOld` - if destiantiong file already exists, the origial will be copied and saved with this string append to the name e.g. `'.bad'`. (Note that if such a file already exists it will not be re-coppied.)      
