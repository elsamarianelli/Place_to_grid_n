function [ spkT ] = genSynData_spkTrain( nFld, sgm, tStp, pkRate, path )
%GENSYNDATA_SPKTRAIN Gen spk train for synthetic blobby field
% Part of project to test false postive rate of gridness based methods.
% This code takes an animals path (real or otherwise) that runs through a
% square environment [size assumed to be 1m^2] then distributes a number
% (nFld) of identical gaussian fields with sigam (sgm) at random. Fields
% can overlap which they do additavly. Then use the path to generate a spk
% train assuming poisson spiking.
%
% ARGS
% nFld          number of fields to distribute
%
% sgm           sigma of fields (they are circular guassians)
%
% tStep         time step in s - needs to be small <1/500s. Note must be an
%               integer divisor into 50hz
%
% pkRate        peak firing rate of each field. NB actuall rate may be
%               if fields coincide
%
% path          xy points [nPts x 2] sampled at 50hz. Normalised to be
%               within 0 to 100 range (i.e. units are cm).
%
% RETURNS
% spkT          times of spikes (sampled at tStep) so min time increment
%               between adjacent spikes is tStep


% --- 
% House keeping
if rem(0.02, tStp) %Check tStp divides exactly into 50hz
    error('tStp must divide exactly into 50hz (0.02)');
end


% ----
%First generate random centres for the fields - just a flat random
%distribution
rng('shuffle');             %Set the rand generator to random seed

fldCnt      =rand(nFld,2)*100; %xy pair between 0 and 100


% ---
%Second generate firing probability at each point on the path based on the
%sum of pdf relative to each field. Loop over each field then sum
% *** Debug stuff ***
% Use this and second bit of debug to get a ratemap to look at - note is
% scale to give desired peak rate below the next imagesc
% [xx, yy]    =meshgrid(0:100);
% path        =[xx(:), yy(:)];
% ***

frRate      =zeros(size(path,1), nFld);
for nn      =1:nFld
    frRate(:,nn)    =mvnpdf(path, fldCnt(nn,:), [sgm^2, sgm^2]); %normal pdf - takes variance
end
% *** MORE DEBUG ***
% imagesc(reshape(frRate(:,1), size(xx)));
% ***

frRate      =sum(frRate,2); %sum of probs

%The norm pdf generated from mvnpdf is normalised to sum to 1 - hence peak
%is equal to 1/(2*pi*sgm). Set peak to pkRate then multiply by tStp to get
%expected number of spks per bin. Hmm think this needs to be squared
frRate      =frRate.* (2*pi*sgm.^2) .* pkRate .*tStp;

% ---
%Convert expectd spk number to a list of when spikes occur accordign to
%Poisson point process. Note this is done for each time step by generating
%a random number 0 to 1, a spike is considered to be fired if this number
%is <= the expected firing rate. NB. allowed states is either a spike
%occurs or no spike occurs, multiple spikes are not allowed but this is
%fine for small tStep (see Dyan & Abott p30)
frRate      =repmat(frRate(:)', [0.02/tStp, 1]); %Convert to tStp sample rate
spkT        =find(rand(size(frRate))<=frRate); %Index with spikes
clear frRate

spkT        =spkT.*tStp; %Spike times


end

