function [ post ] = decode_calcPost_twoPhoton( obsTran, pAgS, varargin)
%DECODE_CALCPOST_TWOPHOTON Given 2p transients & 'ratemap' get posterior probability
%
% Adapted from decode_calcPost with reference to Etter et al 2020 in
% Frontiers. Basic point is to convert the code to work with 2p transients
% which can be treated as a binary event per frame (~frame typicaly at
% 30Hz). So then maths becomes much simpler - for a given cell (ROI) we
% know from exploration data what the probability of observing a transient
% in each spatial bin is p(activity|spatial) aka p(a|s). We want to know
% p(s|a) so that we can decode. For the full Bayes formulism we need p(s)
% and p(a). Because we want a flat prior p(s) is constant for all bins and
% equal to 1/number of bins. p(a) - which can be seen as a normalisation
% term - is just the mean activity over all spatial bins.

% Decodes spatial location based on expected transients and observed 
% transients. Note this assumes that 'transients' have been extract from
% the calcium traces - one way to do this would be to use spike
% deconvolution and asign a 1 or 0 to each frame based on if one is present
% or not
%
% Note1 assumes a flat prior - i.e. not shaped by previous dwell probability
%
% Note2 thid argument (which is optional) allows the normalisation to be
% turned off so the posterior doesn't sum to 1. Be careful with this since
% ratemap bins with 0hz spikes will each yield a probability of 1 if no
% spikes are observed.
%
% Note3 highly likely this code will not work with unvisted (nan) bins -
% was really designed for use in 1D tracks where nans are not encountered
% but will, in prinipal, work in 2D tracks (as long as there are no nans)
%
%
% ARGS          2 or 3 required
% obsTran       [nCells x nFrame] binary (0 or 1) whether transient
%               observed for each cell per frame.
%
% pAgS          [env size 1 x env size 2 x nCells] stacked probability map
%               for all cells being analysed. Specifically probability of 
%               activity given spatial bin or p(a|s). For each bin is the
%               probability of observing a spike per frame. Should be 
%               smoothed with unvisted bins set as nan. Order of cells in 
%               3rd dim should match obsSpks. If 1D track dim2 should have 
%               size=1.
%
%
% normalise     [not required - if not present defaults to true]. Specifiy
%               whether the resulting posterior should be normalised to sum
%               to 1. Either true or false
%
%
% RETURNS
% post          Posterior prob [envSz1 x envSz2 x nFrames] probability
%               animal is in each bin in environment
%
%
% e.g.
% post          =decode_calcPost(mySpikes, myRms, 0.2); %Will normalise
% or
% post          =decode_calcPost(mySpikes, myRms, 0.5, false); %No normalise


% --- HOUSE KEEPING ---
% Need to add a small number to avoid having zeros in the expected number
% of spikes. If expected spikes for a given bin is 0 and we have a single
% spike then the probability of this is 0 which is clearly not useful.
smallV          =eps.^8;

nCell           =size(obsTran,1);
[envSize(1), envSize(2), ~]=size(pAgS);
nFrame          =size(obsTran,2);

if nCell ~= size(pAgS,3)
    error('Number of cells in tBin and rm do not match');
end

if nargin==2 %Proceed as normal with normalisation if not specified
    normalise       =true;
end

if nargin==3 %Normalisation has been specified
    if ~islogical(varargin{1})
        error('3rd variable ''normalise'' must be empty, true or false');
    else
        normalise       =varargin{1};
        clear varargin
    end
end






% --- MAIN ---
% ---
%First stage is to calculate p(a|s), p(s) and p(a)
%Add small number to the p(a|s) otherwise bins with zero transients
%will also yield p(s|a)=0 when a transient is observed.
pAgS              =pAgS + smallV; %[env size 1 x env size 2 x nCells] 

%Second get p(a) probability active for each cell - this is just taken as
%the mean over spatial bins
pA              =nanmean(nanmean(pAgS,1),2); %[1x1xnCell]

%Third use a flat prior over all spatial bins
pS              =ones(size(pAgS,1), size(pAgS,2));
pS              =pS./sum(pS(:)); % [envSize1 x envSize2]

% ---
%Second stage is to calculate the p(s|a) for each cell and spatial bin if
%there was a transient and if there was not a transient
% Now calculate p(s|a) - first calc for trans observed.
pSgA_trans      =bsxfun(@times, pAgS, pS); %[env size 1 x env size 2 x nCells] 
pSgA_trans      =bsxfun(@rdivide, pSgA_trans, pA); %[env size 1 x env size 2 x nCells] 

% And calculate p(s|a) for no trans observed
pSgA_noTrans    =bsxfun(@times, 1-pAgS, pS); %[env size 1 x env size 2 x nCells]
pSgA_noTrans    =bsxfun(@rdivide, pSgA_noTrans, 1-pA); %[env size 1 x env size 2 x nCells] 


% --- 
%Third stage - for each time bin select the appropriate pSgA for each cell
%based on if a transient was observed. A bit clunky but we use a loop over
%frames to do this
 
%Preallocate for speed
tmpPost         =repmat(pSgA_trans, [1,1,1,nFrame]);%tmpPost is 4D [env size 1 x env size 2 x nCells x frames]

for nn          =1:nFrame %Start loop
    currA       =logical(obsTran(:,nn)); %[nCell x1] for current frame 0 or 1 transient from each cell  
    tmpPost(:,:,~currA,nn)  =pSgA_noTrans(:,:,~currA);
end
clear currA
    
    
% ---
%Forth. Take product over all cells to get a p(s|a) for each spatial bin
%i.e. this is a probability density (not normalised) over spatial bins BUT
%to avoid underflow due to very low numbers that will occur due to large
%numbers of ROIs take log first. Note this is a point at which the code
%departs from Etter et al (2020) becuase they use log1p followed by expm1
%partly because they will have some values that are 0 and hence give -inf
%when logged.
%
% NB this still doesn't work - the probabilities are still so low e.g.
% values in post are ~ -2000 before that's taken exp that the resulting
% value is just 0 (i.e. below matlabs threshold). So now revert to the
% Etter method of using log1p and expm1
%
%NB if using the useNoTrans = false then the ROIs with no transients will
%just be left as flat distribution.

tmpPost         =log1p(tmpPost);
post            =sum(tmpPost, 3); 
post            =expm1(post); %Non normalised prob [envSize1 x envSize2 x 1 x nFrames]



post            =reshape(post, [envSize(1), envSize(2), nFrame]);
clear tmpPost

% Finally do we normalise so that the post sums to 1 - this is the default
% position
if normalise
    post        =bsxfun(@rdivide, post, nansum(nansum(post,1),2)); %Normalise 1&2 dim sum to 1
end


end

