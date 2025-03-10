#Analysis_Grid

Code primarily intended for analysing grid cells (e.g. gridness, elipticalness etc). Note many of these functions work on the spatial autocorrelation (SAC) of the ratemap (normally of the smoothed ratemap). A SAC can be generated from a smooth ratemap using the function *xPearson.m* which is in the directory *Analysis_General*. 


## List of contents
Might not be complete - you can help by completing it!

### Directories
**2DFFT_ANL** Folder containing the functions for conducting Fourier analysis of grid cell ratemaps. Basically finds the most significant frequencies composing a spatial firing patter.

### Files
**GRIDELLIPSE_CORRECT** Typically used in conjunction with *gridEllipse_fit.m* - based on the parameters returned by *gridEllipse_fit.m* this function attempts to remove elliptical irregularities from a spatial autocorrelogram (SAC). The regularised SAC can then by assessed for gridness using *multiGridness.m*

**GRIDELLIPSE_FIT** Takes a spatial autocorr (SAC) and fits an ellipse to the six central peaks on the SAC and returns the parameters of that ellipse. Can be used alone to judge how distorted the grid is or in conjunction with *gridEllipse_correct* to remove those distortions prior to calculating gridness with *multiGridness.m*

**MULTIGRIDNESS** Takes a spatial autocorr (SAC) and assess how hexagonally regular it is using two methods. The basic standard version of gridness first used by Sargolini et al 2006: find six central peaks of SAC, extent ROI out beyond these, then assess using correlations how much rot60 symmetry that section has. Also uses the expanding gridness as used in the baby rat grid cell papers. Rather than define a single ROI the region of interest is expanded out and gridness is checked for each radius taking the best value. The absolute values of gridness returned by these methods are NOT directly comparable.



