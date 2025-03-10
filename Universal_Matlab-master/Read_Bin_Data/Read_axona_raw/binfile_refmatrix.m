function [chMat,recDurSamps,RBDLen] = binfile_refmatrix(channels, RBDLen,recDurSamps)

% this creates a reference matrix for extraction of channel's contiguous samples in read_raw_bin_wav.m
% written by M Rutledge, 2/9/2011

remap_channels=[32, 33, 34, 35, 36, 37, 38, 39,0, 1, 2, 3, 4, 5, 6, 7,40, 41, 42, 43, 44, 45, 46, 47,8, 9, 10, 11, 12, 13, 14, 15,48, 49, 50, 51, 52, 53, 54, 55,16, 17, 18, 19, 20, 21, 22, 23,56, 57, 58, 59, 60, 61, 62, 63,24, 25, 26, 27, 28, 29, 30, 31 ];

chMat(1:length(channels),1:recDurSamps)=NaN; 
for channel=channels
    channeled=remap_channels(channel); 
    for i=1:3
        %   16 is packet header. 64 is gap between same channel consecutive
        %   samples within a packet. 216 is packet length
        chMat(channel,i:3:recDurSamps)=16+1+channeled+64*(i-1):216:RBDLen; 
    end   
end