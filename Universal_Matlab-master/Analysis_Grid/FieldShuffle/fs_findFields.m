function [ output_args ] = fs_findFields( binPos, binSpk )
%fs_findFields Part of field shuffle. Find each GC field & surround
% As first step in the field shuffle project take pos and rate data,
% calculate a smooth ratemap and detect the field peaks plus area around
% it. All sections of the ratemap should be allocated to a ratemap.
%
% Step1 is smoothing. Several possible apporaches here:
% i) use adaptive smoothing
% ii) smooth with gaussian kernel or similar
%
% Step2 detect fields and peaks. Again several options:
% i) use watersheding
% ii) detect fields using a threshold, find peak of each field then
% delaunary teselate
%
% DEPENDS
% make_smooth_ratemap CB_Univseral/Analysis_Place_Field

% --- VARS ----------------------------------------------------------------
gSig        =3;   %Gaussian smoothign sigma in bins (bin is about 2cm 
%!!! need to fix bin size in cm!!





% --- MAIN ----------------------------------------------------------------

% First generate a smooth ratemap - use guassian smoothing, turns out to be
% better than adaptive - also make a non-smooth
smthRm          =make_smooth_ratemap(binPos, binSpk, gSig, 'gaus', 'norm');


% Second segregate entire rm using watersheding - note issue with
% watersheding is that it does not allocate the 'ridge' of the watershed to
% a given field, these are left as 0. So use fixWatershed to correct for
% this
smthRm(binPos==0)   =0;     %remove the nans (unvisted bins)
wsAl            =watershed(-smthRm);
wsAl            =fs_fixWatershed(smthRm, wsAl);


%Third within each region of the watershed find the peak bin as well as the
%location of all other bins relative to it (in polar coords)
fldBinLoc       =fs_binMembLoc(wsAl, smthRm); %Returns a cell array


%Now actually reallocate the location of all fields and bins within them -
%returns a mat the same size as smthRm with the mapping of index in smthRm
%to a new index in the shuffled field
staticBin       =find(binPos==0);   %Bins not to move - unvisted ones
shufInds        =fs_shufFields(smthRm, fldBinLoc, staticBin); 


%Test this apply the shuf to binPos and binSpk then smooth 
smthShufRm      =...
    make_smooth_ratemap(binPos(shufInds), binSpk(shufInds), gSig, 'gaus', 'norm');

figure(3)
subplot(2,2,1), imagesc(smthRm);
subplot(2,2,2), imagesc(wsAl);
subplot(2,2,3), imagesc(smthShufRm);
subplot(2,2,4), imagesc(wsAl(shufInds));

sac         =xPearson(smthRm);
[realScale, realGridness]   =autoCorrProps(sac);




%Quick test
[scale, gridness]   =deal(zeros(100,1));
tic
h   =waitbar(0, 'Looping ...');
for nn      =1:1000
    shufInds        =fs_shufFields(smthRm, fldBinLoc, staticBin); 
smthShufRm      =...
    make_smooth_ratemap(binPos(shufInds), binSpk(shufInds), 5, 'boxcar', 'norm');
sac             =xPearson(smthShufRm);
[scale(nn), gridness(nn)] =autoCorrProps(sac);
  waitbar(nn/1000, h)  
end
toc
close(h)

figure(4)
subplot(1,2,1), 
hist(scale), hold on;
plot([realScale, realScale], [0,200], 'r');
hold off
title('Hist of shuffle scale');
xlabel('Scale in bins')
ylabel('Count')

subplot(1,2,2), 
hist(gridness), hold on;
plot([realGridness, realGridness], [0,300], 'r');
hold off
title(['Hist of gridness 95Prctile=' num2str(prctile(gridness, 95))]);
xlabel('Gridness')
ylabel('Count')



end

