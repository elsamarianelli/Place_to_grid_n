function [smthRate, smthSpk, smthPos]=skaggsAdaptiveBinning(pos,spk,alpha)
% Bin spatial data to balance amount of data vs bin size.
% 
% Adaptive binning as used by Skaggs, normally a precursor to calculating
% spatial information.
% 
% TAKES
% pos           unsmoothed binned pos data (in seconds)
% spk           unsmoothed binned spk data
% alpha         smoothing criteria [200]
% sample_rate   rate at which pos sampled [50]

% Orginally from Tom Wills

% Check for empty spk maps %
if sum(sum(spk))==0
    smthPos=pos;    smthPos(pos==0)=nan;
    smthSpk=spk;    smthSpk(pos==0)=nan;
    smthRate=spk;   smthRate(pos==0)=nan;
    return
end
% Pre-assign output %
smthPos=zeros(size(pos));
smthSpk=zeros(size(pos));
% Visited env template: use this to get numbers of visited bins in filter at edge of environemnt %
vis=zeros(size(pos));
vis(pos>0)=1;
% Pre-assign map which records which bins have passed %
smoothedCheck=false(size(pos));
smoothedCheck(pos==0)=true; % Disregard unvisited - mark as already done.

% These parameters depend on place or dir mode %
if size(pos,2)>1
    boundary=0;             % IMFILTER boundary condition
    rBump=0.5;              % Increase radius in 0.5 bin steps.
elseif size(pos,2)==1
    boundary='circular';
    rBump=1;                % Increase radius in 1 bin steps.
end


%%% Run increasing radius iterations %%%
r=1; % Circle radius
while any(any(~smoothedCheck))
    % Check radius isn't getting too big (if >map/2, stop running) %
    if r>max(size(pos))/2
        smthSpk(~smoothedCheck)=nan;
        smthPos(~smoothedCheck)=nan;
        break
    end
    % Construct filter kernel ...
    if size(pos,2)>1
        % Place: Flat disk, where r>=distance to bin centre %
        f=fspecial('disk',r); 
        f(f>=(max(max(f))/3))=1;
        f(f~=1)=0;
    elseif size(pos,2)==1 
        % Direction: boxcar window, r bins from centre symmetrically %
        f=ones(1+(r*2),1);
    end     
    % Filter maps (get N spikes and pos sum within kernel) %
    fSpk=imfilter(spk,f,boundary);
    fPos=imfilter(pos,f,boundary);
    fVis=imfilter(vis,f,boundary);
    % Which bins pass criteria at this radius? %
    warning('off', 'MATLAB:divideByZero');
    binsPassed=alpha./(sqrt(fSpk).*fPos) <= r;
    warning('on', 'MATLAB:divideByZero');
    binsPassed=binsPassed & ~smoothedCheck; % Only get the bins that have passed in this iteration.
    % Assign values to smoothed maps %
    smthSpk(binsPassed)=fSpk(binsPassed)./fVis(binsPassed);
    smthPos(binsPassed)=fPos(binsPassed)./fVis(binsPassed);
    % Record which bins were smoothed this iteration %
    smoothedCheck(binsPassed)=true;
    % Increase circle radius (half-bin steps) %
    r=r+rBump;
end

% Assign Output %
smthRate=smthSpk./smthPos;
warning('on', 'MATLAB:divideByZero');
smthRate(pos==0)=nan;
smthPos(pos==0)=nan;
smthSpk(pos==0)=nan;







