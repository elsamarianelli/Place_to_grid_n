function [gReg, gIReg, gRegShuf, gIRegShuf]=main(gScale, numItt )
%MAIN Test gridness of simulated data bothy gridy and just 'clumpy'
% Writen in response to monkey visual grid paper which appear to have
% spurious grids due to the fact that the null distribution doesn't take
% account of clumpy data.
%
% Code generates regular grids by placing guassian fields in regular array,
% then puts a rats path through this and generates a Poisson spikes train,
% bins to RM and gets gridness. Also does a standard temporal shuffle by
% rotating spikes relative to position. Note the orientation and offset of
% the grid is random on each iteration.
%
% Then does a field shuffle by taking a similar number of fields as would
% appear for regular grids and distributing these around in the window
% randomly. Again generates path etc to get gridness and does the temporal
% shuffle.
%
% Each itteration produces a gridness value for i) reg grid, ii) reg grid
% temp shuffle, iii) ireg grid iv) ireg grid temp shuffle
%
% NB 4 threshs (max due to mem) 12min for 10,000it
%
% ARGS
% gScale - scale of grid to be created in cm [30]
% numItt - number of itterations of code to run
%
% RETURNS
% gReg, gIReg, gRegShuf, gIRegShuf [nIt x 1] vectors of gridness for i)reg
% grid, ii) ireg grid, iii) reg grid temp shuffle, iv) ireg grid temp
% shuffle
%
% e.g.
% [gReg, gIReg, gRegShuf, gIRegShuf]=main(60, 1000 )
%
%
% This function requires several elements to be in the path
% 1) mat file containing the path '1838141210.mat'
% 2) cb_bin_data a wrapper for some of the mTint files
% 3) mTint - required for several functions including bin_pos_data etc



% --- VARS ---------------------------------------------------------------
tStep=1/500; %Time steps in s - needs to be small <1/500
pStep=1/50; %Position sample rate [50Hz]

duration=600; %Lenght of trial in s - must be <= path that is used
peakRate=5; %Peak rate of ideal grid fields-Poisson value may exceed this
baseSigma=0.16; %Sigma for a grid of scale 1


spSamp=0.5; %Spatial sampling freq
smthKern=7; %Smoothing kern for box car smooth of the ratemap

dontGetIdeal=true; %If true saves having to calculate the values for ideal grids

% ---
%Useful precal vars
time2pos=pStep/tStep; %conversion factor between time and pos sampling rate
nTStep=duration/tStep; %Number of time steps
sigma=baseSigma.*gScale; %Sigma for grids of current scale

% --- MAIN ---------------------------------------------------------------
gcp;        %Go parallel

%---- First process path so that it is in form to use
% Load path from a real rat [20 min trial in 1m box] and process
cd('RatPaths');
load('1838141210'); %Creates var pos which is data.pos.*
cd('..');
ppm=str2double(pos.header{strcmp('pixels_per_metre', pos.header),2});
lgthPos=size(pos.xy,1);
if lgthPos/50 < duration, error('Duration must be less than path length'); end
pos.xy=pos.xy-repmat(min(pos.xy), [lgthPos,1]); %get min value to 0,0
pos.xy=pos.xy./ppm*100; %Force values into 0-100 range roughly ppm=100
binSizePix=2;
pos.header{strcmp('window_min_x', pos.header(:,1)),2}='0'; %update header
pos.header{strcmp('window_min_y', pos.header(:,1)),2}='0';
pos.header{strcmp('window_max_x', pos.header(:,1)),2}=num2str(max(pos.xy(:,1)));
pos.header{strcmp('window_max_y', pos.header(:,1)),2}=num2str(max(pos.xy(:,2)));
%Bin pos now
allPos=1:duration/pStep; %Pos points used during trial
posXY = cb_bin_data('dwell_time', 'position', binSizePix, pos, allPos); %binned pos

%---
% Now get path in a format that we can use easily
xPos=interp(pos.xy(:,1), time2pos); %Increase sample rate to match tStep
yPos=interp(pos.xy(:,2), time2pos);
xPos=xPos(1:nTStep); %Truncate path after number of timesteps we want
yPos=yPos(1:nTStep);
xPos=ceil(xPos./spSamp); %now in bins is oversampled space
yPos=ceil(yPos./spSamp);
maxVal=100/spSamp;
xPos(xPos<1)=1; xPos(xPos>maxVal)=maxVal; %Force into bin range - 100 unit env...
yPos(yPos<1)=1; yPos(yPos>maxVal)=maxVal; %... sampled at spSamp
xyInd=sub2ind([maxVal, maxVal], yPos, xPos); %each point on path assigned to bin
clear maxVal xPos yPos


% ---- Now generate field centres for an ideal grid array of unit scale
%Generated grid use scale 1 and orient 0
[xMasterG,yMasterG]=meshgrid(-100:1:100);
yMasterG=reshape(yMasterG.*cosd(30),[],1); % y rows each increase by 0.8660
xMasterG=reshape(xMasterG+repmat(mod([0:0.5:100]',1), [1,201]),[],1); %altenrate x rows offset by 0.5
yMasterG=(yMasterG.*gScale); %Convert to current scale
xMasterG=(xMasterG.*gScale);

%---- Now start loop
[gReg, gIReg, gRegShuf, gIRegShuf]=deal(zeros(numItt,1));
% h=waitbar(0, 'Shuffling ...');
parfor n=1:numItt
    
    %--- Generate base ratemap for ideal grid
    %First get ideal grid with random orient and offset
    rot=rand*60;
    gX=xMasterG + rand*gScale; %Rand offset
    gY=yMasterG + rand*gScale*cosd(30);
    tmpGX=gX.*cosd(rot) - gY.*sind(rot); %Rotation
    gY=gX.*sind(rot) + gY.*cosd(rot);
    gX=tmpGX;
    
    %Final values that are in 200 unit square arena - will ultimatly only
    %use a 100 unit area
    inInd=gX>=0 & gX<=200 & gY>=0 & gY<=200;
    gX=gX(inInd);
    gY=gY(inInd);
    nCents=length(gX); %number of fields in 200 x 200 space
    if ~dontGetIdeal
        [rMap]=getRMap(gX, gY, sigma, peakRate, spSamp); %Regular ratemap
        rMap=rMap.*tStep; %now number of spikes in each tstep
        regExpSpkSeq=rMap(xyInd); %Apply path to map
    end
    
    
    %--- Generate base ratemap for shuffled grid
    nSpk=0;
    while nSpk<20 %Require at lest 20 spikes in ireg so all peaks not off edge
        gX=rand(nCents,1).*200;%Rand x and y pos for peaks - same number ...
        gY=rand(nCents,1).*200; %... as in reg distriubte over 300square
        [rMap]=getRMap(gX, gY, sigma, peakRate, spSamp); %Irreg ratemap
        rMap=rMap.*tStep;
        iregExpSpkSeq=rMap(xyInd);
        
        %Convert expectd spk number to a list of when spikes occur accordign to
        %Poisson point process. Note this is done for each time step by generating
        %a random number 0 to 1, a spike is considered to be fired if this number
        %is <= the expected firing rate. NB. allowed states is either a spike
        %occurs or no spike occurs, multiple spikes are not allowed but this is
        %fine for small tStep (see Dyan & Abott p30)
        %First for reg
        if ~dontGetIdeal
            regSpk=regExpSpkSeq>=rand(size(regExpSpkSeq)); %Logical spikes
            regSpk=find(regSpk); %Index of spikes sampled at tStep
            regSpk=ceil(regSpk/time2pos); %index of spikes sampled at pStep [50hz]
        end
        %Second for ireg
        iregSpk=iregExpSpkSeq>=rand(size(iregExpSpkSeq)); %Logical spikes
        iregSpk=find(iregSpk); %Index of spikes sampled at tStep
        iregSpk=ceil(iregSpk/time2pos);
        nSpk=length(iregSpk); %how many spikes?
    end
    
    %At the same generate temporal shuffled spike train for each of these cells
    %by rotating spikes relative to path note data is alread in ind into pos
    %sampled at 50hz - min shuff +-20s (i.e 1000bins)
    if ~dontGetIdeal
        regSpkShuf=mod(round(regSpk+1000+ rand*(duration*50-2000)), duration*50);
        regSpkShuf(regSpkShuf==0)=1;
    end
    iregSpkShuf=mod(round(iregSpk+1000+ rand*(duration*50-2000)), duration*50);
    iregSpkShuf(iregSpkShuf==0)=1;
    
    %---
    %Now bin and get a ratemap for all four cells
    if ~dontGetIdeal
        spkXY = cb_bin_data('dwell_time', 'position', binSizePix, pos, regSpk); %binned spks
        [~, ~, rmReg] = smooth_field_plot(posXY, spkXY, smthKern,'boxcar'); %Make smooth rm
        spkXY = cb_bin_data('dwell_time', 'position', binSizePix, pos, regSpkShuf); %binned spks
        [~, ~, rmRegShuf] = smooth_field_plot(posXY, spkXY, smthKern,'boxcar'); %Make smooth rm
    end
    
    spkXY = cb_bin_data('dwell_time', 'position', binSizePix, pos, iregSpk); %binned spks
    [~, ~, rmIReg] = smooth_field_plot(posXY, spkXY, smthKern,'boxcar'); %Make smooth rm
    spkXY = cb_bin_data('dwell_time', 'position', binSizePix, pos, iregSpkShuf); %binned spks
    [~, ~, rmRegIShuf] = smooth_field_plot(posXY, spkXY, smthKern,'boxcar'); %Make smooth rm
    
    %Finally get gridness for each
    if ~dontGetIdeal
        sac=xPearson(rmReg);
        [ gReg(n)]=getGridness(sac);
        sac=xPearson(rmRegShuf);
        [gRegShuf(n)]=getGridness(sac);
    end
    sac=xPearson(rmIReg);
    [ gIReg(n)]=getGridness(sac);
%     [~, gIReg(n)]=autoCorrProps(sac);
    sac=xPearson(rmRegIShuf);
    [ gIRegShuf(n)]=getGridness(sac);
%     [~, gIRegShuf(n)]=autoCorrProps(sac);
    
    %     waitbar(n/numItt,h);
end
% close(h);

% ------SUBFUNCTIONS
function [rMap]=getRMap(xCents, yCents, sigma, peakRate, spSamp)
% Given the x and y pos of each gaussian combines them then finds firing
% rate in each sample bin. Note takes pos samples in 0 to 200 range but
% returns a rMap that covers spatial units 50 to 150 in 2d (but sampled at
% rate of spSamp)
% NB note singles to save memory - with large number of fields mem is issue
%
% ARGS
% xCents & yCents [nFields,1] coord of fields centres
% sigma - width of gaussian fields
% peak rate [x] in hz
% spSamp - spatial sample size for generating Poission process
nCents=length(xCents);
[xx,yy]=meshgrid(single(0:spSamp:200));
nElDim=length(xx);
xx=repmat(xx,[1,1,nCents]);
yy=repmat(yy,[1,1,nCents]);
xCents=repmat(shiftdim(single(xCents(:)),-2), [nElDim, nElDim,1]);
yCents=repmat(shiftdim(single(yCents(:)),-2), [nElDim, nElDim,1]);

%Then turn into firing fields
rMap=peakRate.*exp(-  (((xx-xCents).^2)./ (2.*sigma.^2) + ... %GlC field is gaus
    ((yy-yCents).^2)./ (2.*sigma.^2) ));
rMap=sum(rMap,3); %sum and collapse
range=((50+spSamp)/spSamp):1:(150/spSamp);
rMap=double(rMap(range, range)); %Just take block from centre


