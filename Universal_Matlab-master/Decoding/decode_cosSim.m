function [ cosSim, rmCosSim, decBin ] = decode_cosSim( obsSpk, rm)
%DECODE_COSIM Cosine similarity between expected and observed spk vectors
% Frovided with a stacked population ratemap and observations of the same
% spikes - calculats the cosine similarity (0 to 1 - 0 being no
% similarity) between the spikes observed in each time bin and all firing
% vectors expected from each spatial bin. Can be used as a means of
% decoding i.e. pick the bin wit the best match. Advantage is that cosSim
% is not affected by the magnitude of the vectors (over all firing rate)
% e.g. if the cosim of x,y =z then 2x,y = z
%
% NB because code is agnostic of overall rate changes time bin is not a
% required variable.
%
% NB2 bins that have all zero activity in (either observed or expected)
% will generate a nan when compared to.
%
% NB3 IMPORTANT. obsSpk should be observed for time bins that are long
% enough not to introduce too much noise. Would think 1s minimum is a good
% place to start (or smooth the values first).
%
%
% ARGS          3 or 4 required
% obsSpks       [nCells x nTimeBin] number of observed spikes for each cell
%               per time bin. Note will work for non integer values so can 
%               be smoothed.
%
% rm            [env size 1 x env size 2 x nCells] stacked ratemaps for all
%               cells being analysed. Should be smoothed with unvisted
%               bins set as nan. Order of cells in 3rd dim should match
%               obsSpks. If 1D track dim2 should have size=1.
%
%
% RETURNS
% cosSim        [1xnTimeBin] For each time bin the cosine similarity
%               between the observed population vectir and the expected
%               population vector from the spatial bin to which it most
%               closely matches.
%
% rmCosSim      3D [envSize1 x envSize2 x nTimeBins] For each time bing the
%               cosine similarity to the expected vectors from each bin in
%               the ratemap.
%
% decBin        [1xnTimeBin] The spatial bin to which we 'decode' in each
%               timebin i.e. for each time bin the index of the bin in rm
%               with the highest cosSim.
%               
%
%
% e.g.
% cosSim          =decode_calcPost(mySpikes, myRms,); %
%

% --- HOUSE KEEPING ---
[envSize(1), envSize(2), ~]=size(rm);
nTBin           =size(obsSpk,2);


% --- MAIN ---
%Cosime similarity is essentially the dot product of normalised vectors 
%(i.e. unit vectors).

% - First process expected vectors i.e. pop ratemap for each bin - turn
% into a unit vector
%rm is [env size 1 x env size 2 x nCells]
%Calculate the norm (Euclid length) for each spatial bin then divide the
%vector in each spatial bin by that to convert to a unit vector (length=1)
%note no longer have to use bsxfun to do the expansion for element multip.
rmNorm          =sqrt(sum(rm.^2,3)); %[env size 1 x env size 2 x 1]
rmUnit          =rm ./rmNorm;%[env size 1 x env size 2 x nCells]
clear rmNorm


% - Second process observed spikes in a similar way so that everytime bin
% is a unit vector. ObsSpk is [nCell x nTimeBin]
obsSpkNorm      =sqrt(sum(obsSpk.^2)); %[1 x nTimeBin]
obsSpkUnit      =obsSpk ./obsSpkNorm; %[nCell x nTimeBin]
obsSpkUnit          =shiftdim(obsSpkUnit, -2); %3d [1 x 1 x nCells x nTimeBins]
clear obsSpkNorm


% - Third caculate the cosine similarity which is now just the sum of
% product of these normalised vectors between all points in environment and
% the observed rate.
wrking        =rmUnit .* obsSpkUnit; %[envSize1 x envSiz2 x nCells x nTimeBins]
rmCosSim      =sum(wrking,3); %[envSize1 x envSize2 x 1 x nTimeBins]
rmCosSim      =reshape(rmCosSim, [envSize(1), envSize(2), nTBin]); %3D [envSize1 x envSize2 x nTimeBins]

% - Fourth get for each time bin the best cosSim and the bin it came from
[cosSim, decBin] =max(reshape(rmCosSim, [prod(envSize), nTBin]));



end



