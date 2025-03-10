%Script to interpret the output of genSynData_masterLoop
%Data produces by that function are:
%grid           [nScales x nExemplars, types of gridness] gridness of each
%gridFrmShuff   [nScales x nExemplars x nShufTyes x nShuffs of each]
%pVal           [nScales x nExemplars x nShufTypes] 95th prctile of shuff

scales          =20:5:70;

% scales          =[40,80]
nScales         =length(scales);

% ----
% -- TOP ROW - SHOW GRIDNESS FOR EACH OF FOUR METHODS AND ON EACH PLOT
% ILLUSTRATE THE 95PRCTILE

%First plot mean normal gridness vs scale
subplot(2,4,1)
plot(scales, nanmean(grid(:,:,1),2),'k'); %Normal gridness
hold on
plot(scales, nanmean(pVal(:,:,1,1),2), 'r-');
plot(scales, nanmean(pVal(:,:,1,2),2), 'g-');
xlabel('Equivalent grid scale');
ylabel('Grid score');
axis square
title('Norm g: mean with 95% for temporal (red) and field (green)');
hold off

%Second plot Brandon normal gridness vs scale
subplot(2,4,2)
plot(scales, nanmean(grid(:,:,2),2),'k'); %Normal gridness
hold on
plot(scales, nanmean(pVal(:,:,2,1),2), 'r-');
plot(scales, nanmean(pVal(:,:,2,2),2), 'g-');
xlabel('Equivalent grid scale');
ylabel('Grid score');
axis square
title('Brandon g: mean with 95% for temporal (red) and field (green)');
hold off

%Third plot mean expanding gridness vs scale
subplot(2,4,3)
plot(scales, nanmean(grid(:,:,3),2),'k'); %Normal gridness
hold on
plot(scales, nanmean(pVal(:,:,3,1),2), 'r-');
plot(scales, nanmean(pVal(:,:,3,2),2), 'g-');
xlabel('Equivalent grid scale');
ylabel('Grid score');
axis square
title('Expanding g: mean with 95% for temporal (red) and field (green)');
hold off

%Forth plot Brandon expanding gridness vs scale
subplot(2,4,4)
plot(scales, nanmean(grid(:,:,4),2),'k'); %Normal gridness
hold on
plot(scales, nanmean(pVal(:,:,4,1),2), 'r-');
plot(scales, nanmean(pVal(:,:,4,2),2), 'g-');
xlabel('Equivalent grid scale');
ylabel('Grid score');
axis square
title('Brandon expanding g: mean with 95% for temporal (red) and field (green)');
hold off


% ----
% -- BOTTOM ROW NOW SHOW THE PROPORTION OF BLOBBY CELLS CLASSIFED AS GRIDS
% AT 95PRC
%Fifth plot norm grid false postive rates
subplot(2,4,5)
plot(scales, sum(grid(:,:,1) > pVal(:,:,1,1),2)./size(pVal,2), 'r');
hold on
plot(scales, sum(grid(:,:,1) > pVal(:,:,1,2),2)./size(pVal,2), 'g');
plot(scales, ones(nScales,1)*0.05, 'k-');
xlabel('Equivalent grid scale');
ylabel('Proportion blobby cells classified as grids');
ylim([0, 0.2]);
title('Norm grid: Red vs temporal shuffle, green vs field shuffle')
axis square
hold off

%Sixth plot Brandon grid false postive rates
subplot(2,4,6)
plot(scales, sum(grid(:,:,2) > pVal(:,:,2,1),2)./size(pVal,2), 'r');
hold on
plot(scales, sum(grid(:,:,2) > pVal(:,:,2,2),2)./size(pVal,2), 'g');
plot(scales, ones(nScales,1)*0.05, 'k-');
xlabel('Equivalent grid scale');
ylabel('Proportion blobby cells classified as grids');
ylim([0, 0.2]);
title('Brandon grid: Red vs temporal shuffle, green vs field shuffle')
axis square
hold off

%Seventh Expanding grid false postive rate
subplot(2,4,7)
plot(scales, sum(grid(:,:,3) > pVal(:,:,3,1),2)./size(pVal,2), 'r');
hold on
plot(scales, sum(grid(:,:,3) > pVal(:,:,3,2),2)./size(pVal,2), 'g');
plot(scales, ones(nScales,1)*0.05, 'k-');
xlabel('Equivalent grid scale');
ylabel('Proportion blobby cells classified as grids');
ylim([0, 0.2]);
title('Expanding grid: Red vs temporal shuffle, green vs field shuffle')
axis square
hold off

%Eighty plot expanding Brandon grid false postive rates
subplot(2,4,8)
plot(scales, sum(grid(:,:,4) > pVal(:,:,4,1),2)./size(pVal,2), 'r');
hold on
plot(scales, sum(grid(:,:,4) > pVal(:,:,4,2),2)./size(pVal,2), 'g');
plot(scales, ones(nScales,1)*0.05, 'k-');
xlabel('Equivalent grid scale');
ylabel('Proportion blobby cells classified as grids');
ylim([0, 0.2]);
title('Brandon expanding grid: Red vs temporal shuffle, green vs field shuffle')
axis square
hold off
