function swap_list = led_swap_filter(led_pos, led_pix)
% Checks for instances of two leds swapping or big one replacing little one
% when the big one gets obscured.
% Input xy posiiton of each led and
% and number of pixels in each. Big light is light number 1
% format: led_pos(1:n_pos, 1:num_cols, x-y), npix(1:n_pos, 1:num_cols)
%
% This original version of this function can be found at this end of this
% file inside a comment block.  The current version is more effecient at
% handling nans and thus faster.
%
thresh = 5;

%get nan mean and nan std of led_pix
led_pix_tmp = led_pix;
are_nans = isnan(led_pix_tmp);
led_pix_tmp(are_nans) = 0;
n_ok = [size(led_pix,1) size(led_pix,1)]- sum(are_nans,1);
mean_npix = sum(led_pix_tmp,1)./n_ok;

denom = max(n_ok-1,1);
denom(n_ok==0)=nan;
led_pix_tmp = bsxfun(@minus,led_pix_tmp,mean_npix);
led_pix_tmp(are_nans) = 0;
std_npix = sqrt(sum(led_pix_tmp.^2,1)./denom);
%------------------------------------


% Check if big light closer to small light at t-1 than to big light at t-1
% and small light either not found or closer to big light at t-1 than small light at t-1

% Use one of the two following blocks of code - the first calculates a city
% block metric and the second euclidsian distance. For most applications the
% latter is correct.

% Calculate city block metric
% dist12 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 2, 2));
% dist11 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 1, 2));
% dist21 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 1, 2));
% dist22 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 2, 2));

%Calculate eucldian
led_pos_tmp = led_pos(1:size(led_pix,1),:,:);
led_pos_tmp(isnan(led_pos_tmp)) = 0;


dist12 = sqrt(sum(((led_pos_tmp(2:end,1,:)-led_pos_tmp(1:end-1,2,:)).^2),3));
dist11 = sqrt(sum(((led_pos_tmp(2:end,1,:)-led_pos_tmp(1:end-1,1,:)).^2),3));
dist21 = sqrt(sum(((led_pos_tmp(2:end,2,:)-led_pos_tmp(1:end-1,1,:)).^2),3));
dist22 = sqrt(sum(((led_pos_tmp(2:end,2,:)-led_pos_tmp(1:end-1,2,:)).^2),3));

pos = 2:size(led_pix,1);
switched = (dist12 < dist11-thresh) & ( isnan(led_pos(pos, 2, 1)) |(dist21 < dist22-thresh) );

% Check if size of big light has shrunk to be closer to that of small light (as Z score)
z11 = (mean_npix(1) - led_pix(2:end, 1))/std_npix(1);
z12 = (led_pix(2:end, 1) - mean_npix(2))/std_npix(2);
shrunk = z11 > z12;
swap_list = find( switched & shrunk ) + 1;

% ------------------ The original version of the function follows ---------------------------

%{

function swap_list = led_swap_filter(led_pos, led_pix)
% Checks for instances of two leds swapping or big one replacing little one
% when the big one gets obscured.
% Input xy posiiton of each led and
% and number of pixels in each. Big light is light number 1
% format: led_pos(1:n_pos, 1:num_cols, x-y), npix(1:n_pos, 1:num_cols)

thresh = 5;
%ok_pos = find(all(~isnan(led_pix),2));

mean_npix = nanmean(led_pix);
std_npix = nanstd(led_pix);

% Check if big light closer to small light at t-1 than to big light at t-1
% and small light either not found or closer to big light at t-1 than small light at t-1
pos = 2:size(led_pix,1);

% Use one of the two following blocks of code - the first calculates a city
% block metric and the second euclidian distance. For most applications the
% latter is correct.

% Calculate city block metric
% dist12 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 2, 2));
% dist11 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 1, 2));
% dist21 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 1, 2));
% dist22 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 2, 2));

%Calculate eucldian
dist12 = sqrt(nansum(((squeeze(led_pos(pos,1,:))-squeeze(led_pos(pos-1,2,:))).^2),2));
dist11 = sqrt(nansum(((squeeze(led_pos(pos,1,:))-squeeze(led_pos(pos-1,1,:))).^2),2));
dist21 = sqrt(nansum(((squeeze(led_pos(pos,2,:))-squeeze(led_pos(pos-1,1,:))).^2),2));
dist22 = sqrt(nansum(((squeeze(led_pos(pos,2,:))-squeeze(led_pos(pos-1,2,:))).^2),2));

switched = (dist12 < dist11-thresh) & ( isnan(led_pos(pos, 2, 1)) |(dist21 < dist22-thresh) );

% Check if size of big light has shrunk to be closer to that of small light (as Z score)
z11 = (mean_npix(1) - led_pix(pos, 1))/std_npix(1);
z12 = (led_pix(pos, 1) - mean_npix(2))/std_npix(2);
shrunk = z11 > z12;
swap_list = find( switched & shrunk ) + 1;

% ----------------------------------------------------------------------------------------------------------------------

%}
