function [ isGrid, grid, pVal ] = identify_grid_cell( pos, spkPos )
%IDENTIFY_GRID_CELL Identify grid by comparing gridness to shuffle
%
% ARGS
% pos           The the pos section of the data structure loaded by
%               read_DACQ e.g. data.pos which contains pos.xy, pos.header
%               etc
%
% spkPos        Pos samples at which condiate cells spikes were recorded

% RETURNS
% isGrid        [4x2] logical indicating if the cell is a grid at 95% CI vs
%               8 ways of doing the shuffle. The first column being
%               standard gridness, eliptical corrected standard gridness,
%               expanding gridness, and eliptical corrected standard
%               gridness. These are all compared against a temporal
%               shuffle. The second column is the same thing but against a
%               field shuffle.
%
% grid          [4x1] the gridness value of the cell based on the four
%               methods described above
%
% pVal          [4x2] Not a p-value, actually the 95% threshold for each of
%               the four methods (rows) and for each of the two shuffles
%               (columns - temporal followed by field).
%
% DEPENDS
% fs_masterAnalysis
% xPearson
% mTint files



% --- VARS
nLpForShuf          =100; %For each rm how many shuffles to do get sig [100]
binSzCm             =2;    %Bin size in cm [2]


% --- MAIN

%First work out bin sizes
ppm             =str2num(pos.header{strcmp('pixels_per_metre', pos.header(:,1)),2});
binSize         =(ppm/100)*binSzCm; %Pixels per 2cm


%Second do the main analysis - gets the gridness for cell then does 100
%shuffles.
[ grid, pVal, ~, ~  ] ...
                =fs_masterAnalysis( data.pos, spkPos, binSize, nLpForShuf);

            
%Third determine if this is a grid agains the 95% of the shuffles -
isGrid          =repmat(grid, [1,2]) > pVal;
end

