function [ grid, pVal, gridFrmShuff ] = genSynData_masterLoop(  )
%GENSYNDATA_MASTERLOOP Gen syntetic grid data, test if qualifies as grid
% Master function for the field shuffle project. A series of loops through
% which we generate 'blobby' data i.e. spatial multi peaked fields that
% aren't periodic (drop 2d gaussians randomly). Generate poison spike
% trains, bin and get ratmap then test vs own shuffle to see if the cells
% qualify as grids. Aim is to try and bring the false postive rate down
% with a better shuffle. Loop is over different grid scales.
%
% To get path I have taken real rat paths and normalised to be between 0
% and 100 then cut into 1min chunks. (saved in posChunks.mat). Idea is to
% generate 20minute (or however long) paths by randomly recombining these.
% posChunks is a 3D mat of size [nPosSampPerChunk x 2 x nChunk] contains xy
% pairs normally for 1minute chunks of data
%
% I basically calculate four versions of gridness: standard gridness which
% I've always used, expanding gridness (described in Langston et al 2010)
% and to each of these I apply the Brandon correction for elipticity to
% give another two versions. Each of these four are compared against values
% from a shuffled data set (two shuffled for each grid - temporal and
% spatial). Hence for each of the four types of gridness I have two null
% distruibutions
%
% RETURNS
% grid      3D mat [nScales x nLoop x 4] for each set of synthetic data the
%           gridness measured by the four different gridness measures
%           (looped over scale and a number of iterations in each scale)
%
% pVal      4D mat[nScales x nLoop x 4 x 2] similar to above but also with
%           the 95th confidence interval for each of the four types of
%           gridness calculated for temporal then spatial shuffle. NOT A
%           PVALUE
%
% gridFrmShuff  4D mat [nScales x nLoop x 4 x nShuf] for all the synethetic
%           data the shuffles based on their ratemaps. Can be used to get
%           pVal

%Set the middle - kk - loop to be parfor
fprintf('Remember to set the middle loop to be parfor [line 172] ...\n');

% ---
%VARS
nLpPerScale         =1000; %For each scale how many random blobby cells to test [1000]
nLpForShuf          =100; %For each rm how many shuffles to do get sig [100]
gScl2Test           =30:10:80; %Grid scales to test in cm [20:10:80]

%Debug settings
% nLpPerScale         =5; %For each scale how many random blobby cells to test [100]
% nLpForShuf          =50; %For each rm how many shuffles to do get sig [1000]
% gScl2Test           =[20,40]; %Grid scales to test in cm

%Next determine the field width - modelled by 2D gaussian. This value will
%be used to get sigma (the std) of the gaussian. Note variance is sigma^2.
%95% of gaussian is in 2xsigma from centre. If we define 95% as field edge
%and we want the edge of neighbouring fields to be a full filed width
%appart (i.e. 4*sigma) then sigma = gScale/8
sg2Gsc              =0.1250; %Ratio of sigma to grid scale sgm = gScl * sg2Gsc

%How should number of fields scale with the gScl if grid was cartesian
%would be ^2 of number of cycles in one dim. But if we assume grid is
%aligned with x axis then have 100/gScal cycles in that axis but more in
%the y axis since the y distance to the next axis parallel to x is 0.8660
%wavelenghts - tand(60)/2. Hence get more cycles in the y dim 2/tand(60).
%Total number of filds is numCyclesPerDim^2 * 2/tand(60)
nFld2Gsc            =2/tand(60); %nFlds = round((100/gScl).^2 * nFld2Gsc)

tStp                =1/500; %Must be small <=1/500 and exaclty divide into 1/50

trlLgth             =20; %Number of 1min path chunks to use
pkRate              =10; %Peak firing rate

binSize             =2; %Bin size in cm

%NB smoothing options only for gridness calc - field seg is always based on
%gaus
smthKernSz          =5; %smooth kern size - bins for boxcar and sigma for Guass [5 or 3]
smthKernTyp         ='boxcar'; %'boxcar' or 'gaus'

minTmpRot           =30; %Minimum temporal rotation in s

% ---
% HOUSE KEEPING
gcp             %Start pool if not started
load posChunks; %Loads the mat with poschunks in


% --- Start main
%Preallocate to store gridness, pvalue from dist and entire distribution
nGScl2Test          =length(gScl2Test);
%NL grid is the gridness of the data a 3D mat with [nScales, nExamples, nTypesG]
% Currently that 3rd dim has length 4 since there are four types of Gridness
% Standard, Standard with Brandon Correction, Expanding and Expanding with 
% Brandon Correction
grid                =zeros(nGScl2Test, nLpPerScale,4); %Grid of data
pVal                =zeros(nGScl2Test, nLpPerScale,4,2); %Store p val for each type of gridness
gridFrmShuff        =zeros(nGScl2Test, nLpPerScale,4,2,nLpForShuf);


% ----
%Some precalcs that are used in loops below - don't want to redo everytime
trlLgthPosSamp      =trlLgth * 60 *50; %Trial length in post sampl
%For temporal shuffle - calculate the allowed shift in pos sample
rotRng          =[minTmpRot, trlLgth*60 - minTmpRot]*50; %min & max rot in possamples

%cb_bin_data needs the full pos struct - put some dummy values in
pos.header      ={'sample_rate', '50'; ... %Required for bin_pos
    'window_max_x', '100';...
    'window_min_x', '0';...
    'window_max_y', '100';...
    'window_min_y', '0'};
pos.speed       =20; %Ditto
pos.dir         =2;
pos.xy          =[];

% ---
% Loop over grid scales
h                   =waitbar(0, 'Looping over grid scales ...');
for nn              =1:nGScl2Test
    gScl            =gScl2Test(nn); %Current grid scale
    
    % Loop over multiple synthetic cells in each scale
    for kk              =1:nLpPerScale
        
        %For each cell - first get a path by picking random chunks and
        %concatanating
        chnk2Use        =randperm(size(posChunks,3));
        chnk2Use        =chnk2Use(1:trlLgth); %Random path
        path            =posChunks(:,:,chnk2Use); %still 3D
        %Now convert to a nPosSamp x 2 mat
        path            =reshape(permute(path, [2, 1, 3]), 2, []);
        path            =path';
        
        nFld            =round((100./gScl).^2 * nFld2Gsc); %Number fields
        spkT            =genSynData_spkTrain( nFld, gScl*sg2Gsc, tStp, pkRate, path );
        spkPos          =ceil(spkT .*50); %Pos samp for each spike
        
        
        % ---
        %Now have spike train and and path - bin these to get a ratemap
        %then calculate gridness as required.
        %First bin
        pos.xy          =path;
        
        binXY           =cb_bin_data('dwell_time', 'position', binSize, pos, 1:length(path));
        binSpk          =cb_bin_data('spikes', 'position', binSize, pos, spkPos);
        
        
        rm              =make_smooth_ratemap(binXY, binSpk, smthKernSz, smthKernTyp, 'norm');
        rmGaus          =make_smooth_ratemap(binXY, binSpk, 3, 'gaus', 'norm');
        sac             =xPearson(rm);
        [grid(nn,kk,1), grid(nn,kk,3)]...
            =fast_get_gridness( sac ); %Gridness of the data stdG & expG
        
        %Brandon correction - see subfunction - does eliptical correct
        [grid(nn,kk,2), grid(nn,kk,4)]...
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
            rotBinSpk       =cb_bin_data('dwell_time', 'position', binSize, pos, rotSpkPos);
            rm              =make_smooth_ratemap(binXY, rotBinSpk, smthKernSz, smthKernTyp, 'norm');
            sac             =xPearson(rm);
            [tmpTSG, tmpTEG]=fast_get_gridness( sac ); %Gridness of the data
            [tmpTSBG, tmpTSEBG]=sf_getBrandonG(sac); %Brandon g on temporal
            
            
            
            %2) Do the new fields shuffle
            shufInds        =fs_shufFields(wsRm, fldBinLoc, staticBin);
            rm              =...
                make_smooth_ratemap(binXY(shufInds), binSpk(shufInds), smthKernSz, smthKernTyp, 'norm');
            sac             =xPearson(rm);
            [tmpSSG, tmpSEG]=fast_get_gridness( sac ); %Gridness of the data
            [tmpSSBG, tmpSSEBG]=sf_getBrandonG(sac); %Brandon g
            
            %The 8 types of gridness are (ignoring tmp):
            %TSG temporal shuffle standard G
            %TEG temporal shuffle expanding G
            %TSBG temporal shuffle Brandon correction G
            %TSEBG temporal shuffle expanding with Brandon correct G
            %S*** is then the same sequence but for Spatial shuffle
            gridFrmShuff(nn,kk,:,:,hh)     =[tmpTSG, tmpTEG, tmpTSBG, tmpTSEBG; ...
                                            tmpSSG, tmpSEG, tmpSSBG, tmpSSEBG]';
            
        end %Finish loop over shuffle
        clear tmp*
        
        %Get 95th conf int based on shuf
        pVal(nn,kk,1,:)   =prctile(squeeze(gridFrmShuff(nn,kk,1,:,:)),95,2);
        pVal(nn,kk,2,:)   =prctile(squeeze(gridFrmShuff(nn,kk,2,:,:)),95,2);
        pVal(nn,kk,3,:)   =prctile(squeeze(gridFrmShuff(nn,kk,3,:,:)),95,2);
        pVal(nn,kk,4,:)   =prctile(squeeze(gridFrmShuff(nn,kk,4,:,:)),95,2);
        
    end
    waitbar(nn/nGScl2Test,h);
end
close(h)

end

%--- SUBFUNCTION
function        [bGStd, bGExp]      =sf_getBrandonG(sac)
% Get brandon gridness - do his correction then get gridness
        [ rSac1, rSac2 ] =brandonRegEcc( sac);
        [tmp1Std, tmp1Exp]=fast_get_gridness(rSac1);
        [tmp2Std, tmp2Exp]=fast_get_gridness(rSac2);
        bGStd       =max(tmp1Std, tmp2Std); %take best
        bGExp       =max(tmp1Exp, tmp2Exp); %take best
end

