function circDiv     =dir_kldiv(dirRM)
%Returns a kl divergence for directional firing: measure of directionality.
% Calculates Kullback-Leibler diveregence between a directional ratemap 
% and a circular distribution of the same mean. Hence is a measure of how
% directional a cells firing is. Note should really be run on a smoothed 
% polar ratemap - since the measure has no concept of the continuity of the 
% bins in a polar distribution lack of smoothing/undersampling will yield
% high KL scores.
%
% Tom Wills found that kl divergence correlates very highly with bits/spike 
% (Skaggs info). The measures are probably almost identical. He also found 
% a threshould of >0.25 was a reasonable cut off (cells above this being 
% cosidered to be directional). Though in papers such as Doeller et al
% (2010) we used >0.15.
%
% NB. calls kldiv - function from MathWorks file exchange - must be in path
%
% TAKES
% dirRM     Binned and smoothed directional ratemap
%
% RETURNS
% circDiv   [scalar] KL divergence from circular distribution


%VARS
%Can't define KL divergence for values of 0 (will return nan as log2(0) is -inf. So if
%there are zeros add a very small number to the zero values.
smallIncr       =0.00001;


dirRM(dirRM==0) =smallIncr;

normSmoothDir   =dirRM./sum(dirRM(:)); %Normalise to a probability
nDirBins        =length(dirRM);
compCircle      =ones(1, nDirBins)/nDirBins;
circDiv         =kldiv((1:nDirBins)', normSmoothDir(:), compCircle(:));