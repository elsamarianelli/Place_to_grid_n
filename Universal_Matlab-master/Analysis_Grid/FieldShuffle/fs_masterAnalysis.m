function [ grid, pVal, gridFrmShuff, gridScl, allRm] = fs_masterAnalysis( pos, spkPos, binSize, nLpForShuf)
%FS_MASTERANALYSIS Encapsulates all the shuffles & gridness
% One function that contains or calls all code required to build ratemaps,
% get gridness and do the shuffle. The idea being this single function can
% be applied to synethic or real data - hence ensures exactly the same
% analysis is used and saves on coding.
%
% So for a single cell (real or otherwise) this generates the ratemap and
% sac, gets gridness (four versions) then shuffles some number of times
% specified by nLpForShuf (both temporal and field), gets the gridness of
% each of those and returns the ranking of the real gridness against that.
%
% ARGS
% pos       pos part of data strucutre returned by read_DACQ ie. pos =
%           data.pos
%
% spkPos    pos sample that each spike is assigned to (i.e. index into
%           pos.xy.
%
% binSize   spatial bin size in pix, should generally be set to yield 2cm
%           bin (hence depends on ppm)
%
% nLpForShuf Number of shuffles to do - both temporal and field
%
%
%
% RETURNS
% grid      Gridness of the real cell (4 versions)
%
% pVal      The 95th confidence interval for each of the 4 types of
%           gridness for each of the two types of shuffle
%
% grdFrmShff Stores the gridness of each shuffle iteration for reference
%
% gridScl   [scalar] Grid scale in bins of the actual SAC
%
% allRm     [optional] A 3 cell cell array containing. cell{1} the original
%           ratemap, cell{2} another cell array containing all the time
%           shuffled ratemaps, cell{3} another cell array containing all
%           the spatial shuffled ratemaps. If this output isn't specified
%           then these ratemaps are not saved and the code will be faster.


% --- VARS ----------------------------------------------------------------
gcp;
%NB smoothing options only for gridness calc - field seg is always based on
%gaus
smthKernSz          =1.5; %smooth kern size - bins for boxcar and sigma for Guass [5 or 3]
smthKernTyp         ='gaus'; %'boxcar' or 'gaus'

minTmpRot           =30; %Minimum temporal rotation in s

if nargout==5 %We need to return the ratemaps
    keepRm          =true;
else
    keepRm          =false;
end


% --- MAIN ----------------------------------------------------------------
%Pre allocate
grid                =zeros(4,1); %Grid of data
pVal                =zeros(4,2); %Store p val for each type of gridness
gridFrmShuff        =zeros(4,2,nLpForShuf);
if keepRm
    timeShuffRm     =cell(nLpForShuf,1);
    spaceShuffRm    =cell(nLpForShuf,1);
end



binXY           =bin_data('dwell_time', 'position', binSize, pos, 1:size(pos.xy,1));
binSpk          =bin_data('spikes', 'position', binSize, pos, spkPos);

rm              =make_smooth_ratemap(binXY, binSpk, smthKernSz, smthKernTyp, 'norm');
if keepRm
    allRm{1}    =rm;
end
rmGaus          =make_smooth_ratemap(binXY, binSpk, 3, 'gaus', 'norm');
sac             =xPearson(rm);
[grid(1), grid(3), gridScl]...
    =fs_fastGridness( sac ); %Gridness of the data stdG & expG

%Brandon correction - see subfunction - does eliptical correct
[grid(2), grid(4)]...
    =sf_getBrandonG(sac); %Again stdG & expG both with correction




% ---
%Some preparation for the field shuffle - working on the ratemap
%find the fields and segment them
wsRm            =rmGaus;
wsRm(binXY==0) =0;     %remove the nans (unvisted bins)
wsAl            =watershed(-wsRm);
wsAl            =fs_fixWatershed(wsRm, wsAl);
fldBinLoc       =fs_binMembLoc(wsAl, wsRm); %Returns a cell array
staticBin       =find(binXY==0);   %Bins not to move - unvisted ones


%Some precalcs that are used in loops below - don't want to redo everytime
trlLgthPosSamp      =size(pos.xy,1); %Lenght of trial in pos samples
trlLgth             =trlLgthPosSamp/(50*60); %Lenght of trial in mins
%For temporal shuffle - calculate the allowed shift in pos sample
rotRng          =[minTmpRot, trlLgth*60 - minTmpRot]*50; %min & max rot in possamples

% ---
%Finally do the shuffles
%NL whould be parfor
parfor hh          =1:nLpForShuf
    
    % ----
    %1) Do the standard temporal shuffle - rotate spikes by some
    %number of pos points - faster than doing on the actual spike
    %times - though does keep spikes together in a possample
    rotAmnt         =ceil(rand(1) * (rotRng(2) - rotRng(1)) + rotRng(1)); %in pos samples
    rotSpkPos       =mod(spkPos + rotAmnt, trlLgthPosSamp);
    rotSpkPos(rotSpkPos==0)=nLpForShuf;
    rotBinSpk       =bin_data('spikes', 'position', binSize, pos, rotSpkPos);
    rm              =make_smooth_ratemap(binXY, rotBinSpk, smthKernSz, smthKernTyp, 'norm');
    if keepRm
        timeShuffRm{hh}     =rm;
    end
    sac             =xPearson(rm);
    [tmpTSG, tmpTEG]=fs_fastGridness( sac ); %Gridness of the data
    [tmpTSBG, tmpTSEBG]=sf_getBrandonG(sac); %Brandon g on temporal
    
    
    
    %2) Do the new fields shuffle
    shufInds        =fs_shufFields(wsRm, fldBinLoc, staticBin);
    rm              =...
        make_smooth_ratemap(binXY(shufInds), binSpk(shufInds), smthKernSz, smthKernTyp, 'norm');
    if keepRm
        spaceShuffRm{hh}     =rm;
    end
    sac             =xPearson(rm);
    [tmpSSG, tmpSEG]=fs_fastGridness( sac ); %Gridness of the data
    [tmpSSBG, tmpSSEBG]=sf_getBrandonG(sac); %Brandon g
    
    %The 8 types of gridness are (ignoring tmp):
    %TSG temporal shuffle standard G
    %TEG temporal shuffle expanding G
    %TSBG temporal shuffle Brandon correction G
    %TSEBG temporal shuffle expanding with Brandon correct G
    %S*** is then the same sequence but for Spatial shuffle
    gridFrmShuff(:,:,hh)     =[tmpTSG, tmpTSBG, tmpTEG, tmpTSEBG; ...
        tmpSSG, tmpSSBG, tmpSEG, tmpSSEBG]';
    
end %Finish loop over shuffle
clear tmp*

%Get 95th conf int based on shuf
pVal(1,:)       =prctile(squeeze(gridFrmShuff(1,:,:)),95,2);
pVal(2,:)       =prctile(squeeze(gridFrmShuff(2,:,:)),95,2);
pVal(3,:)       =prctile(squeeze(gridFrmShuff(3,:,:)),95,2);
pVal(4,:)       =prctile(squeeze(gridFrmShuff(4,:,:)),95,2);

if keepRm
    allRm{2}        =timeShuffRm;
    allRm{3}        =spaceShuffRm;
end

end


%--- SUBFUNCTION
function        [bGStd, bGExp]      =sf_getBrandonG(sac)
% Get brandon gridness - do his correction then get gridness
[ rSac1, rSac2 ] =brandonRegEcc( sac);
[tmp1Std, tmp1Exp]=fs_fastGridness(rSac1);
[tmp2Std, tmp2Exp]=fs_fastGridness(rSac2);
bGStd       =max(tmp1Std, tmp2Std); %take best
bGExp       =max(tmp1Exp, tmp2Exp); %take best
end