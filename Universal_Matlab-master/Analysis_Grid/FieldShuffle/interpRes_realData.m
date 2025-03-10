%allData the cell array should be loaded
%The basic point is that I have calculated 8 gridness thresholds - based on
%these have 95% and compare with real values
%
% Data is stored in cell array 'allData' , columns are as follows:
% 1) filename and path
% 2) tetrode
% 3) cell number
% 4) grid scale in cm
% 5) all four gridness scores
% 6) 95% prctile for all 8 gridness scores (temp shuff first then spatial)
% for each: stdG, expG, stdG Brandon, expG Brandon.
% 7) Full list of shuffled g


%Loop over every row of allData (i.e. each cell)
nCell       =size(allData,1);

isGrid      =zeros(4,2,nCell); %Logical to record if each cell is a grid
oldIsGrid   =zeros(nCell,1);
for nn      =1:nCell
   %Loop down and check if each is above it's own threshold 
    isGrid(:,:,nn)      =repmat(allData{nn,5}, [1,2])> ...
                        allData{nn,6};
    oldIsGrid(nn)       =allData{nn,5}(1)>=0.3;
   
end

%Finally get the proportion of grids based on each method
isGrid      =reshape(isGrid, [8,nCell]);
propGrid    =mean(isGrid,2);
propOldGrid =mean(oldIsGrid);










fprintf(['/nTotal number of cells: ' num2str(nCell)]);
fprintf('/nGrid seq is temporal first, stdG, expG, stdG Brandon, expG Brandon the spatial');
propGrid
propOldGrid

