function [ p, d, smth_p, smth_d ] = wrapper_pxd( spkInd, data, bin_size )
%WRAPPER_PXD Calls bins data, calls pxd
% Simple wrapper for pxd code - abstracts out the binning up of direction
% and pos data, also does some basic plotting.
%
% ARGS
% spkInd        pos index corresponding to each spike (i.e. as found in
%               data.tetrode(x).pos_sample).
%
% data          full mtint data strucutre (e.g. carries data.pos,
%               data.tetrode etc.). Important bit here is actually data.pos
%
% bin_size      [1 x 3] The bin sizes to use for dir, posY, then posX 
%               ([6,16,16]).
               




% ---
% 1) Bin the data into the 3d mat required by pxd. Specifically pxd needs
% to have [dir,posY,posX] mat.
spkDirPos       =[data.pos.dir(spkInd), data.pos.xy(spkInd,[2,1])];
allDirPos       =[data.pos.dir(:), data.pos.xy(:, [2,1])];

%Determine the data ranges including minimum values
dataRng         =[360, max(data.pos.xy(:,2)), max(data.pos.xy(:,2)); ...
                0, min(data.pos.xy(:,2)), min(data.pos.xy(:,1))];

spkBin         =histnd(spkDirPos, dataRng, bin_size);
allBin         =histnd(allDirPos, dataRng, bin_size)./50; %in seconds


% ---
% 2) Now call the pxd code
[p, d]          =pxd(spkBin, allBin);


% ---
%3) Smooth and plot - note uses incorrect smoothing as I'm smoothing across
%the unvisted bins. Really these should be ignored and set to nan.

%First do ratemap
unvisited       =squeeze(sum(allBin,1))==0; %unvisited bins
smthKern        =ones(5)./5; %boxcar
denom           =filter2(smthKern, ~unvisited);

tmp_p           =p;
tmp_p(unvisited)=0;
smth_p          =filter2(smthKern,tmp_p);
smth_p          =smth_p./denom;
smth_p(unvisited)=nan;

%Then polar for which there are no unvisited dir
smthKern        =ones(3,1)./3;
smth_d          =imfilter(d, smthKern, 'circular');

%Plot
subplot(1,2,1)
imagesc(smth_p);
axis equal tight
subplot(1,2,2)
theta           =linspace(0,2*pi,(360/bin_size(1))+1)';
polar(theta(1:end-1),smth_d);

end

