function [ meanAng ] = cb_circMean( angs, varargin )
%CB_CIRCMEAN Circular mean - to replace unreliable circ_mean
% Mean over circular values (in rads).
%
% ARGS - take 1, 2 or 3 [second two optional]
% angs - mat of angles in rads - mean is over first non-single dim
% dim - dimension to work along, if not supplied uses first non singelton
% w - weighting of angles - if not supplied assumed to be equal
%
% e.g.
% [meanAng] = cb_circMean(someAngs);
% [meanAng] = cb_circMean(someAngs, [], dim2average);
% [meanAng] = cb_circMean(someAngs, someWeights, dim2average);

% --- HOUSE KEEPING ------------------------------------------------------
switch nargin
    case 1 %Just supply angs
        dim=find(size(angs)>1,1); %1st non sing
        w=ones(size(angs));
        
    case 2 %Supply angs and dim
        dim=varargin{1};
        w=ones(size(angs));
        
    case 3 %Supply angs, weights and dim to work on
        if isempty(varargin{1}) %Check if [] given
            dim=find(size(angs)>1,1); %1st non sing
        else
            dim=varargin{1};
        end
        w=varargin{2};
        
    otherwise %Error message
        error('Supply one or two arguments');
end
clear varargin


% --- ERROR CHECKS -------------------------------------------------------
if size(angs)~=size(w)
    error('Size of weight and ang args must be the same');
end


% --- MAIN ---------------------------------------------------------------

%Make sure nans are mutual
w(isnan(angs))=nan;
angs(isnan(w))=nan;

sumW=nansum(w, dim);

meanAng=atan2(nansum(sin(angs).*w, dim)./sumW, nansum(cos(angs).*w, dim)./sumW);
meanAng=mod(meanAng,2*pi); %Values in 0 to 2pi range

end

