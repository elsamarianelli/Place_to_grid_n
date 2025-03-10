function shufInds       =fs_shufFields( smthRm, fldBinLoc, staticBin)
%FS_SHUFFIELDS Shuffles field locs & attempts to keep relative loc of bins
% Functinon that actually does buisness. Randomly reallocates the position
% of the field centres. Then works through the fields starting in each with
% the bin that is closest to the centre, replaces this bin relative to the
% new field centre location (NB deals with the fields in a random
% sequence). Then moves to the next closest bin for each field again
% attempting to preserve the closeness. NB2. Because several bins will have
% the same distance to the field centres works round in a CCW direction but
% for each field starts at a random orientation.
%
% ARGS
% smthRm        Smooth ratemap (only used to get the size of the ratemap)
%
% fldBinLoc     Cell array [nField x 5] returned from fs_binMembLoc that
%               specifies the location of each field centre and the bins
%               relative to that.
%
% staticBin     List of indicies not to move (typically unvisted bins) i.e
%               find(posBin==0)
%
% RETURNS
% shufInds      A 2d mat same size as smthRm with all ind reallocated to
%               new locations following rules described above

% --- SOME BASIC HOUSE KEEPING ------------------------------------------
nField          =size(fldBinLoc,1);
bin2mv          =1:numel(smthRm);
bin2mv(staticBin)=[]; %Bins available to use
nBin            =length(bin2mv);
shufInds        =zeros(size(smthRm)); %preallocate the returned arg


%reorder the rows of fieldBinLoc to be a random sequence - this will be the
%order we deal with the fields - don't want there to be a fixed sequence
fldBinLoc(randperm(nField),:) ...
                =fldBinLoc;  %reorder fieldBinLoc following fieldOrd


% --- FIRST REALLOCATE THE FIELD CENTRES - THE EASY STEP -----------------
fldCentInd       =bin2mv(randperm(nBin));
fldCentInd       =fldCentInd(1:nField);

%Now start to populate shufInd - first with new field centre locations then
%add in the unvisted bins that will not be moved (staticInds)
shufInds(fldCentInd)=[fldBinLoc{:,2}];
shufInds(staticBin)=staticBin;



% --- NOW DECIDE ON THE ORDER THE EACH BIN IN EACH FIELD WILL BE DELT WITH-
%Bins from each field are dealth with taking the nearest first but apply a
%rule to deal with tied field: take bin with smallest angle to peak first
%(all values should be 0 to 2pi) so this would mean we move CCW from
%horizontal. But want to add a degree of randomness so for each field add a
%fixed random orientation. Essentially this is equivalent to rotating each
%field before reallocating it.

rndOrient       =rand(nField,1) .*(2*pi); %Random orient to add to each fld

for nn          =1:nField %Have to loop to work through cell array
    %First remove bins listed against each field that are static and should
    %not be moved
    [~, ind2rmv]        =intersect(fldBinLoc{nn,3}, staticBin);
    fldBinLoc{nn,3}(ind2rmv)=[];
    fldBinLoc{nn,4}(ind2rmv)=[];
    fldBinLoc{nn,5}(ind2rmv)=[];
    fldBinLoc{nn,6}     =fldBinLoc{nn,6} - length(ind2rmv); %update number of bins
    
    %Add random orient
    fldBinLoc{nn,4}     =mod(fldBinLoc{nn,4}+rndOrient(nn), (2*pi));

    %Now define the sequence to take bins in based first on dist to peak
    %then orientation (starting with closest to 0)
    if fldBinLoc{nn,6}==0 %Field is only one bin - the peak 
        continue
    end
    [~, tmpInd]         =sortrows(cell2mat(fldBinLoc(nn,4:5)),[2,1]);
    tmpIndSeq           =fldBinLoc{nn,3}; %Current sequence of ind in fld
    fldBinLoc{nn,3}     =tmpIndSeq(tmpInd);%Put them back in the new seq.
end
clear tmp*


% ---  ORDER BINS ACCORDING TO CLOSENESS TO NEW FIELD CENTRES ------------
%For every bin centre find the distance and angule to all bins in the
%ratemap - will use this as the basis for reallocation
[tmpI, tmpJ]    =ind2sub(size(smthRm), 1:numel(smthRm)); %subscript for all bins row


%So get subscript locations of bin centres and subtract from subscript locs
%of all bins
[tmpFldCI, tmpFldCJ]=ind2sub(size(smthRm), fldCentInd);
tmpI            =bsxfun(@minus, tmpI, tmpFldCI'); %[nField x nBin]
tmpJ            =bsxfun(@minus, tmpJ, tmpFldCJ');
[tmpTh, tmpR]   =cart2pol(tmpJ, tmpI); %[nField x nbin]

%As before add a random orient to the bins related to each field -
%basically equivalent of starting to place data into bins at a random
%offset from each field centre. Need to do this to avoid systematic biases.
rndOrient       =rand(nField,1) .*(2*pi); %Random orient to add to each fld
tmpTh            =mod(bsxfun(@plus, tmpTh, rndOrient), (2*pi));
clear rndOrient

%Finally have to loop over the fields to sort the bins in the order in
%which to deal with them. Aim is to create a 2D mat [nField x nBin] which
%is for each field centre (row) the index of all bins in the rm in the
%order of closeness to the field centre. The final step is to set the bins
%that have already been taken (by field centres) to 0 (so they don't get
%reused)
binOrdByFld     =zeros(size(tmpTh));
for nn          =1:nField
   [~,tmpInd]       =sortrows([tmpR(nn,:)', tmpTh(nn,:)']);
   [~, tmptmp]      =intersect(tmpInd, fldCentInd); %Field centres are taken
   tmpInd(tmptmp)   =0;
   [~, tmptmp]      =intersect(tmpInd, staticBin); %Static bins don't move.
   tmpInd(tmptmp)   =0;
   binOrdByFld(nn,:)=tmpInd(:);
end
clear tmp*

% --- NOW THE REASON WE CAME HERE - REALLOCATE ALL BINS TO NEW LOCATIONS --
%This has to be an iterative loop, can't be vector or parfor. Basic method
%is iterate round the rows in fldBinLoc, the 3rd column contains cells with
%the bin indices in each field arranged in order of closeness to the peak;
%select these in turn. Each one is then place in the bin specified by the
%row in binOrdByFld matching the current cell. The bin to use is then the
%first non-zero ind in that row. Set it to zero after along with it's
%occurance in all other rows.

%First have to pull out the sequence of ind that will be dealt with (i.e.
%take one from first field, move to next field etc). Have to loop over
%field to do this
nBinPerFld      =cell2mat(fldBinLoc(:,6)); %nBins in each field
[tmpInd, tmpId] =deal(ones(nField, max(nBinPerFld)) .*nan); %Preal
for nn          =1:nField
    tmpInd(nn,1:nBinPerFld(nn))=fldBinLoc{nn,3};
    tmpId(nn,1:nBinPerFld(nn))=nn; %
end
%Extract just the non-nans
binIndSeq       =tmpInd(~isnan(tmpInd)); %Sequence of bins to shuff ...
binFldIdSeq     =tmpId(~isnan(tmpInd));  % and the fields they belong to
clear tmp*

%So finally do the shuffle
for nn      =1:length(binIndSeq) %number of bins to shuf (nBins - fld cent)

    %Decide which field we're currently working with
    curFldInd     =binFldIdSeq(nn); %Current bin belongs to this field
    
    %Now determine which is the next unallocated bin close to this field
    %centre
    nxtBin      =binOrdByFld(curFldInd,find(binOrdByFld(curFldInd,:),1));
    
    %Update shufInd to indicate the new location of this bin
    shufInds(nxtBin)=binIndSeq(nn);
    
    %Finally update binOrdByFld to indicate that this bin is taken and
    %should not be used again. This is done for all rows (i.e. for all
    %fields)
    binOrdByFld(binOrdByFld==nxtBin)=0;
end
 
end