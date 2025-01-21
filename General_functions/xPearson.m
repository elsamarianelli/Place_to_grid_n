function [xCoef, nBins, lags]=xPearson(seq1, varargin) 

% Pearson correlation between two sequences at all lags.
% Conducts a Pearson correlation between two vectors or matricies at all possible lags
% (i.e. 1d or 2d). So broadly similar to xcorr but returns Pearson correlation coefficient
% not some unspecified coefficient. Can deal with nans in a ratemap, these are just
% ignored.
% NB. Tried to use fastfilter2 which used FFT to do fast conv but for very
% small numbers (e.g. 0) would often produces values like 10x10-12 which
% lead to issues later (e.g. corrs >1)
%
% TAKES
% seq1 - a vector or matrix, if only one value is supplied does an autocorrelation
% seq2 - supply a second vector or matrix to do a cross correlation.
% filterType - value passed to filter2 ('full' default, 'same' or 'valid')
%
% RETURNS
% xCoef - pearson corr for all lags [similar dimensions to seq1 & seq2 e.g. for autocorr
%           of a row mat [1xn] returns [1x2n-1]]
%
% nBins - for each lag the number of bins used to calculate the result (if nans are
%           present these are not counted in the total)
%
% lags      [same size as xCoef] for each postion the offset in the input
%           sequences used to calcluate that data point. NB central point
%           is always 0
%
% EXAMPLE
% [xCoef, nBins, lags]=xPearson(sind(0:719))
%
% NB
% For a 1d example, using sind(2*(0:359)), I've tested the results for 0 to 4 lag against
% corrcoef and get the same numbers.


if nargin==1
    seq2=seq1;
    filterType='full';
elseif nargin==2
        seq2=varargin{1};
         filterType='full';
elseif nargin==3
            seq2=varargin{1};
         filterType=varargin{2};
else
    error('Too many arguments');
end


%Deal with nans in seq1 or seq2, these should be ignored (i.e. don't contribut to nBins
%and shouldn't be proliferated by filter2)
onesSeq1    =ones(size(seq1));
onesSeq1(isnan(seq1))=0; %Set nans to be 0 in ones mat
onesSeq2    =ones(size(seq2));
onesSeq2(isnan(seq2))=0; %Set nans to be 0 in ones mat

seq1(isnan(seq1))=0; %Set nans to be 0 in main sequence
seq2(isnan(seq2))=0;

seq1Sq      =seq1.^2;
seq2Sq      =seq2.^2;


seq1xseq2   =filter2(seq1, seq2, filterType); %Returns mat of size(seq1)+size(seq2)-1 [equivalent to xcorr(seq1, seq2)
sumSeq1     =filter2(seq1, onesSeq2, filterType);
sumSeq2     =filter2(onesSeq1, seq2, filterType);
sumSeq1Sq   =filter2(seq1Sq, onesSeq2, filterType);
sumSeq2Sq   =filter2(onesSeq1, seq2Sq, filterType);
nBins       =filter2(onesSeq1, onesSeq2, filterType);
nBinsSq     =nBins.^2;



warning('off'); %Need this as we get quite a few divide by zeros
stdSeq1     =(sumSeq1Sq./nBins - (sumSeq1.^2)./nBinsSq).^0.5;
stdSeq2     =(sumSeq2Sq./nBins - (sumSeq2.^2)./nBinsSq).^0.5;

covar       =seq1xseq2./nBins - (sumSeq1.*sumSeq2)./nBinsSq;

xCoef       =covar./(stdSeq1.*stdSeq2);
xCoef       =real(xCoef); %Remove any complex which sometimes appear

warning('on');

%Finally caculate lags if it is required
if nargout==3
    %Lags required
    
    lgthCoef =size(xCoef,2);
    cntrPt   =(lgthCoef/2)+0.5;                                             %centre bin of xCoef
    lags     =(1:lgthCoef)-cntrPt;                                          %the lags
end

end
