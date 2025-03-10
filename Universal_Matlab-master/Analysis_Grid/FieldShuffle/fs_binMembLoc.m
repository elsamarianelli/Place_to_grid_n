function [ fldBinLoc ] = fs_binMembLoc( wsAl, rm)
%FS_BINMEMBLOC Find peak in each field & asign dist&orient of each bin
% As part of the field shuffle have already found the fields using
% watersheding. Now for each field find the peak bin and for every other
% bin in the field define, in polar, coords the location relative to the
% peak. Will try to preserve these during the field shuffle.
%
% ARGS
% wsAl          Assignment of all bins to a field as a mask. Dim is same as
%               as ratemap with interger values.
%
% rm            The smooth ratemap all this is based on
%
% RETURNS
% fieldBinLoc   A cell array [nField, 4] So rows are different fields. 
%               Then columns are as follows: 1) Field mask number, 2) Index
%               of peak  3) index of each other bin 4)theta to peak for 
%               each bin other than peak bin 4) distance to peak for each 
%               bin

szRm            =size(rm); 

% First find the peaks in each field - this is the easy step.
unqFieldMask    =unique(wsAl(:)); %Different labels used
nField          =length(unqFieldMask);%number of fields

% Second do the work with a loop - cycle through each field and get the
% peak plus polar coord to all other bins in the field. Save these in a
% cell array [nField, 4] So rows are different fields. Then columns are as
% follows: 1) Field mask number, 2) Index of peak  3) index of each other
% bin 4)theta to peak for each bin other than peak bin 5) distance to peak 
% for each bin
fldBinLoc       =cell(nField, 4); %Preallocate

for nn          =1:nField
   mask             =wsAl==unqFieldMask(nn);
   tmp              =rm;
   tmp(~mask)       =nan; %Just the current field outside is nan
   [~,maxInd]       =max(tmp(:)); %Index of max val in field
   otherInd         =find(mask); %All ind in field
   otherInd         =otherInd(otherInd~=maxInd); %Ind of other bins in field
   
   %Now calc offset between otherInd and maxInd and convert to polar coord
   %with origin along x-axis theta increasing anti-clock
   [maxI, maxJ]     =ind2sub(szRm, maxInd);
   [otherI, otherJ] =ind2sub(szRm, otherInd);
   offSetI          =otherI-maxI;
   offSetJ          =otherJ-maxJ;
   [offSetTh, offSetR]=cart2pol(offSetJ, -offSetI);
   
   %Put it all in the cell array
   fldBinLoc{nn,1}  =unqFieldMask(nn);  %Field mask number
   fldBinLoc{nn,2}  =maxInd;            %Index of peak bin of field
   fldBinLoc{nn,3}  =otherInd;          %Index of each other bin
   fldBinLoc{nn,4}  =offSetTh;          %Angle to other bins in field
   fldBinLoc{nn,5}  =offSetR;           %Dist to other bins in field
   fldBinLoc{nn,6}  =length(otherInd);  %Number of other bins
end



end

