function [ scale2test, falsePos, thresh95 ] = wraper_scale_range(  )
%WRAPER_SCALE_RANGE Wraper for main - tests lots of dif grid scales
% Runs main function with different scales to get false positve rate of
% temporal shuffle vs field shuffle
%
% NB for small scales with large numbers of fields this will cause out of
% memory if running on many workers (suggest limit to 4 workers min scale
% 20)
% --- VARS ---------------------------------------------------------------
scale2test=[10:5:80]; %try 10:5:80 [10:5:80] - can just get away with 4 workers
nIt=10000; %Num itt for each step [10000]


% --- MAIN CODE ----------------------------------------------------------
nScale=length(scale2test);

h=waitbar(0, 'Testing different scales...');
falsePos=zeros(nScale,1);
for n=1:nScale
    %Do shuffle using main
[~, gIReg, gRegShuf, gIRegShuf]=main(scale2test(n), nIt );
 thresh95=prctile(gIRegShuf,95); %95th perctile from temporal shuff [use ideal or not]
falsePos(n)=sum(gIReg>=thresh95)./nIt; %Percentage exceeding thresh
    waitbar(n/nScale,h);
end
close(h);
   

plot(scale2test, falsePos); %Draw
end

