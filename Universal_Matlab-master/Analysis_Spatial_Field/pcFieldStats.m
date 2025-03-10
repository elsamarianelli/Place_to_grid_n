function [ cellOn, fieldArea, meanRate, maxFieldArea ] = pcFieldStats( rmS, minField, minRate )
%PCFIELDSTATS Basic metrics about place cell fields
% Takes a smooth ratemap, finds valid fields and returns metric relating to these,
% specifically: field size in bins, mean in field rate. Field is defined as continguous
% area above 2x mean rate and occupying >= minField bins
%
% ARGS
% rmS - smoothed ratemap, unvisted bins as nans
% minField - minimum number of bins for field to be valid
% minRate - minimum mean infield rate to be considered active

% RETURNS
% cellOn 1 or 0 - indicates where cell meets criteria for being active (has field and
%       rate)
% fieldArea - area of valid fields in bins
% meanRate - mean in field rate (so mean of bins in the field)
%
% E.G.
% [ cellOff, fieldArea, meanRate ] = pcFieldStats( rmS, 15, 1 )


%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Code is setup as a loop with a single itteration - this is because there are several
%criteria by which a cell might be counted as 'off', encountering one of the criteria
%resuts in a break from the loop

%%% VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldThresh=2 * nanmean(rmS(:)); %Fields are bins > mean rate



for n=1 %Loop with single it
    %Trivial case: no firing at all, cell off
    if nansum(rmS)==0,
        [cellOn, fieldArea, meanRate, maxFieldArea]=deal(0);
        break, 
    end
    
    %Define valid fields - above field threshold and larger than min size
    labelRm=rmS>fieldThresh; %logic - bins in field
    labelRm=bwlabel(labelRm,4);% connected bins
    stats=regionprops(labelRm,  'area');
    areas=[stats.Area]; %area of each field
    clear stats
    validFields= find(areas>= minField);
    if isempty(validFields),%No valid fields
     [cellOn, fieldArea, meanRate, maxFieldArea]=deal(0);
        break, 
    end 
    
    %There are some valid fields so get combined area and rate
    fieldArea=sum(areas(validFields)); %Combined area of valid fields
    meanRate=mean(rmS(ismember(labelRm, validFields))); %Mean infield rate - valid fields
    maxFieldArea=max(areas(validFields));
    if meanRate<=minRate, cellOn=0; break, end%Firing rate too low to be valid
    
    cellOn=1;
end


end

