function [n_jumpy, led_pos] = led_speed_filter(led_pos, max_pix_per_sample, led)
% Filters out short runs of data caused by tracker picking up an incorrect distant point.
% Resets led_pos_x to NaN
%
% There are two versions of this function. Switch between the two version by
% commenting/uncommenting line 18 and 19. They should both give the same
% output, but the _fast version is significantly more difficult to read. 

%Assume first non-nan is okay
ok_pos = find( ~isnan(led_pos(:,led,1)) );

if length(ok_pos) < 2
    warning(' < 2 tracked points for led %d\n', led);
end
mpps_sqd = max_pix_per_sample^2;

%###### Uncomment one of the following two lines: ########################
%[n_jumpy, led_pos] = led_speed_filter_slow(led_pos, mpps_sqd, led, ok_pos);
[n_jumpy, led_pos] = led_speed_filter_fast(led_pos, mpps_sqd, led, ok_pos);
%#########################################################################
end


function [n_jumpy, led_pos] = led_speed_filter_slow(led_pos, mpps_sqd, led, ok_pos)
%This is the original, simple, but unoptimised version.

n_jumpy = 0;
prev_pos = ok_pos(1);

for p = 2:length(ok_pos)
    pos = ok_pos(p);
    % Get speed of shift from prev_pos in pixels per sample (squared)
    pix_per_sample_sqd = ((led_pos(pos,led,1)-led_pos(prev_pos,led,1))^2+...
        (led_pos(pos,led,2)-led_pos(prev_pos,led,2))^2)/(pos-prev_pos)^2;

    if pix_per_sample_sqd > mpps_sqd
        led_pos(pos, led, :) = NaN;
        n_jumpy = n_jumpy+1; 
    else
        prev_pos = pos;
    end
end

end

function [n_jumpy, led_pos] = led_speed_filter_fast(led_pos, mpps_sqd, led, ok_pos)
%This version is around 5 times faster (depending on data), but it is more complicated.
%
%It has been optimised by noting the following:
%Most of the data will be ok, i.e the difference between successive positions 
%will be smaller than the tollerance. In the cases where it is not small
%enough, skipping one or two samples is almost always enough.  In the small
%number of cases where this is still not enough, we can revert to the original slow 
%algorithm.  But in a dataset with 100,000 pos we might expect about 4 or 5 
%instances where we have to use the slow algorithm.  And each of these
%"problem" points will only be a few tens of pos long.  Thus the worst part 
%of the loop is reduced from having to run 100,000 times down to only about
%a hundred.
%Of the data which is not considered a "problem", we still need a way of
%selecting which points have been "jumped over" and which have been "landed
%on".  It seems very difficult to vectorise this calculation entirely
%(although if we simplified the problem to only permit jumping over one pos
%it would be relatively easy to solve.).  The algorithm used here finds
%blocks of pos which are easy to solve and "summarises" them so that the
%loop, when we get to it, can jump through them quickly.  This takes care
%of about 99.99% of pos.  About 0.01% is actively, but quickly dealt with 
%in the loop, and a tiny number of points are dealt with actively, but
%slowly using the original simple algorithm.
%
%There are more comments in the code below.
%
% DM

led_pos_sqeezed = squeeze(led_pos(ok_pos,led,:));

%----Is a jump of 1 small enough?-------------------------------
diff_t = diff(ok_pos);
diff_xy = diff(led_pos_sqeezed);
jmp = [int8(sum(diff_xy.^2,2)./diff_t.^2 <= mpps_sqd) ; int8(0)];

%----If not, is a jump of 2 small enough?-----------------------
not_jmp = find(~jmp(1:end-2));
diff_t = ok_pos(not_jmp) - ok_pos(not_jmp+2);
diff_xy = led_pos_sqeezed(not_jmp,:)-led_pos_sqeezed(not_jmp + 2,:);
jmp(not_jmp) = int8(2)*int8( sum(diff_xy.^2,2)./diff_t.^2 <= mpps_sqd );

%-----If not, is a jump of 3 small enough?-----------------------
not_jmp = find(~jmp(1:end-3));
diff_t = ok_pos(not_jmp) - ok_pos(not_jmp+3);
diff_xy = led_pos_sqeezed(not_jmp,:)-led_pos_sqeezed(not_jmp + 3,:);
jmp(not_jmp) = int8(3)*int8( sum(diff_xy.^2,2)./diff_t.^2 <= mpps_sqd );

%-----If not, record a value of 4 which indicates >3 ------------
jmp(~jmp) = 4;

clear diff_t diff_xy not_jmp
%--- At this point jmp tells us how far we need to jump if we find
%ourselves at any given element of ok_pos (except where jmp==4).-----

%-- We now create the array 'nice' which will eventually tell us which
%ok_pos to keep and which to discard.  The first stage is to find elements
%which are obviously going to be "landed on" in the loop.  For example:
% jmp = ...111X.... , we know that X will be landed on (unless jmp==4 comes
% into play, in which case we will make corrections).
L = int32(length(ok_pos));
nice = [false; false; false; 
        ( jmp(3:end-3)==1  & jmp(2:end-4)<3 & jmp(1:end-5)<4) ;
        false ; false]; 
    
%--- Consider the example, jmp = ....112X..., we know that X will not be
% "landed on" in the loop.  So, we create an array 'known' which holds this 
% information, i.e. elements which we KNOW either will or wont be landed on.
% Just by looking at the previous three jumps, we CAN'T always be sure if X
% will be landed on.  e.g. ......231X.., X may or may not be landed on.
known = nice;
not_nice = find(~nice(4:end-2))+3;
known(not_nice) = (mod(jmp(not_nice-1)-2,4)<2 & mod(jmp(not_nice-2),2)==1 & jmp(not_nice-3)<3);
clear not_nice

%--- Now we are going to create jmp2, which is like jmp, but it holds absolute 
%addreses rather than relative jumps. i.e jmp=13121.. becomes jmp2=25466...
%This in itself would not be interesting, what is interesting is we can now
%summarise contiguous blocks of known==true.  We find the last nice==true
%in each block and set all the previous nice==true elements to point to it.  Thus
%however the block of known is entered, it will be possible to jump most of
%the way through.  By this method, we can reduce the number of iterations
%of the loop to a few hundred rather than tens of thousands. -------------

jmp2 = int32(jmp);                        
ind_start_known_block = [int32(find(known(1:end-1)<known(2:end))); L];
ind_end_known_block = [int32(find(known(1:end-1)>known(2:end))); L];
ind_end_known_block = ind_end_known_block - int32(~nice(ind_end_known_block)) ...
                        - int32( ~nice(ind_end_known_block) & ~nice(ind_end_known_block-1)); %we need the last nice element of the block.  If this is not the last element, then it must be the second last or third last.
jmp2_temp(ind_start_known_block,:) = ind_end_known_block;
jmp2_temp(ind_end_known_block) =  -ind_end_known_block + jmp2_temp(ind_end_known_block); %this takes care of blocks of length one
jmp2_temp = cumsum(single(jmp2_temp)); %cumsum is not supported for integer types! What?
jmp2_temp(ind_end_known_block) = jmp2(ind_end_known_block)+ind_end_known_block;

ind_problems = find(jmp2==4);
ind_not_nice = int32(find(~nice));
jmp2(nice) = jmp2_temp(nice);
jmp2(ind_not_nice) = jmp2(ind_not_nice)+ind_not_nice;
jmp2(ind_problems) = -ind_problems;
jmp2(end-2:end) = jmp2(end-2:end).*int32(jmp2(end-2:end)>0); %set to zero any of the last 3 elements which is negative


%---And finally, we can do the dastardly loop, but it should be much faster now.--
p=1;
while p ~= 0
    while p>0
        nice(p) = true;
        p = jmp2(p);
    end
    
    if p<0
        Problem = -p; p = -p+3;
        while p<L && ~nice(p)
            p = p+1;
            nice(p) = sum((led_pos_sqeezed(Problem,:) - led_pos_sqeezed(p,:)).^2) ...
                   ./ (ok_pos(Problem) - ok_pos(p)).^2 <= mpps_sqd;
        end
        p = jmp2(p);
    end
end

%And this is what we came here for....
not_nice = ~nice;
led_pos(ok_pos(not_nice),led,:) = nan;
n_jumpy = sum(not_nice);
end