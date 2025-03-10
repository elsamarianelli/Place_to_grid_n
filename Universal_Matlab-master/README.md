# Universal_Matlab
Common Matlab functions used for analysis (e.g. loading data, making rate maps, identifying grid cells) and modelling. Files are arranged in folders based loosely on what they are for.

**If you are new** are good starting point is the *Read_Bin_Data* folder that contains functions for reading in and processing raw Axona data, it also contains some example data sets that can be used to test other functions on. Probably the most useful first function to know is *read_DACQ.m* which is contained in *Read_Bin_Data* and can be used to load raw Axona data. You can find some example data to load in the *Read_Bin_Data/Example_data* directory.

Most functions have good comments and start with a clear function description indicating what they do and how to use them (the function description can be reached either by typing help [function name] or using the Matlab function editor to look at the first few lines of the function.

## List of Contents ##
This list is not yet complete - you can help by completing it.

**Analysis_Freq** Directory containing code for analysis of oscillatory signals (e.g. theta power in LFP, spiking frequency, phase precession etc).

**Analysis_Grid** A directory containing matlab code primarily intended for analysing the properties of grid cells (e.g. gridness, elipticalness etc).

**Figures_Presentation** Directory containing code that relates to figure generation and presentation in general (e.g. the lab colour map for probability posteriors etc.)

