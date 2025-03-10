function newWsAl    =fs_fixWatershed( rm, wsAl )
%FIXWATERSHED Allocate ridges of watershed fields to closest neighbour
% Watersheding for field detection leaves the ridges unallocated - for some
% purposes these need to be assigned to a field. This code checks each
% unallocated bin and assigns it to the neighbour (8 or 4) with the most
% similar value in the ratemap.
%
% NB tries to do this without loop for speed
%
% ARGS
% rm            smooth ratemap
% wsAl          watershed allocation of fields - unallocated bins are 0
%

% --- Define the kernal i.e. which neighbours to look at
% Define what neighbours to look at - should be either 4 (i.e. share a
% side) or 8 (share a side or a corner).
kern            =zeros(3);
kern([2,4,6,8])   =1; %Set as 4


% --- Now identify the bins to reassign
% First find all unallocated bins
[unAlI, unAlJ]  =find(wsAl==0); %As sub

% Create a 3D stack with one bin in each layer of the stack
unAlSt          =zeros(size(rm,1), size(rm,2), length(unAlI));
unAlInd         =sub2ind(size(unAlSt), unAlI, unAlJ, (1:length(unAlI))');
unAlSt(unAlInd) =1;

% Convolve with the kern to define the neighbourhood of each stack to look
unAlStNeig      =convn(unAlSt, kern, 'same');
% But we don't want the neighbours to include bins that are allocated as 0
% so remove these. To do this first have to create a stack of the field 
% labels, trun to logical and combine with neihbours
wsAlSt          =repmat(wsAl>0,[1,1,length(unAlI)]);
unAlStNeig      =unAlStNeig & wsAlSt;
clear wsAlSt

% --- Then select those bins out of the ratemap and compare them
%Create a similar sized ratemap stack to select from
rmSt            =repmat(rm, [1,1,length(unAlI)]);

%Then select - first the neighbourhood values then the value to compare
%against
rmNeighSt       =rmSt;
rmNeighSt(~unAlStNeig)=nan;

rmValSt         =rmSt;
rmValSt(~unAlSt)=0;
rmValSt         =convn(rmValSt, kern, 'same');

%NL is the difference in rate between each bin that is unallocated and each
%of its neighbours (defined by kern). Just need to assign the bin to the
%lowest of these while ignoring nans
rmDifSt         =abs(rmValSt-rmNeighSt);

% --- Finally reassign the unallocated bins to match the neibhour they
% closest resemble
%Get the index in each plane of the min value - convert to 2d to do this.
%This is, for each unallocated bin, the index of the neighbour bin that has
%the closest rate to it
[~,minInd]      =min(reshape(rmDifSt, [], length(unAlI)));

newWsAl         =wsAl;
newWsAl(sub2ind(size(wsAl), unAlI, unAlJ))=wsAl(minInd(:));

end

