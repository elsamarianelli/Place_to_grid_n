function [xy] = cumrot(xy, dtheta, ixy_other)
% xy is a nx2 array of (x_i,y_i) positions at each time point
% dtheta is a nx1 array of rotation values, (th_i) for each 
% time point. The rotation th_i is taken about the point (x_i,y_i)
% not about the origin. Note that this is why the rotation is
% specified as >d<theta, meaning it's the rotation that occurs
% within that time window.
% dtheta should be expressed in radians.
% 
% If a third input is provided, it should be mx3, it is some "other"
% data to apply the computed rotations to. The second and third
% columns of it are xy data, and the first column is an index into the 
% original xy/dtheta data (i.e. a number from 1...n).
%
% As you can see this function is pretty heavily vectorized, and
% currently takes about 13miliseconds to run on 20mins of pos data.
% That's kind of nice, right?
%
% Each rotation is expressed as a rotation matrix of the form:
%  [a b c;
%   d e f;
%   0 0 1]
% Where:
% a = cos(th)
% b = -sin(th)
% c = x - xcos(th) + ysin(th)
% d = sin(th)
% e = cos(th)
% f = y - xsin(th) - ycos(th)
%
% Such a 3x3 matrix can be applied to the point (x',y') by multiplying
% it by [x' ; y' ; 1]. This will rotate (x',y') around the point (x,y)
% by an angle th.
%
% What we want for each i in 1...n, is a rotation matrix which
% is the product of all the rotation matrices up to that point.
% This is like a cumsum, but for matrix products...this class
% of operation is called a "prefix" operation....google "prefix sum"
% if you want to know more.  There are also more notes at the bottom
% of this file.
%
% Sorry about the name of the function.
%
% DM, May 2015.
%

% zero-pad out to pow2 length, for simplicity, we remove padding at end
% (for 20mins at 50Hz this is 10% padding) 
true_n = length(dtheta);
logn = ceil(log2(true_n));
dtheta(end+1:2^logn) = 0;
xy(end+1:2^logn,:) = 0;

% we need to copy these as we will eventually be modifying xy in place
x_ = xy(:,1);  y_ = xy(:,2);

% get the elements of the individual rotation matrices
cos_theta = cos(dtheta);
sin_theta = sin(dtheta);
a = cos_theta;
b = -sin_theta;
c = x_ - x_.*cos_theta + y_.*sin_theta;
d = sin_theta;
e = cos_theta;
f = y_ - x_.*sin_theta - y_.*cos_theta;

% Now do the prefix-matrix-multiply.
% This takes two loops with a total of 2*log2(n)-1 iterations...

% do first log2(n) steps
for log_s = 0:logn-1
    s = 2^log_s;
    s2 = 2*s;
    
    % pick out the right-hand values. 
    a_r = a(s2:s2:end); b_r = b(s2:s2:end); c_r = c(s2:s2:end);
    d_r = d(s2:s2:end); e_r = e(s2:s2:end); f_r = f(s2:s2:end);
    
    % and the left-hand values
    a_l = a(s:s2:end);  b_l = b(s:s2:end);  c_l = c(s:s2:end);
    d_l = d(s:s2:end);  e_l = e(s:s2:end);  f_l = f(s:s2:end);
    
    %Note that we need to store a copy of the RHVs because we modify them
    %in place during each iteration of the loop but need to know the
    %original values from the start of the iteration.
    
    % store in right hand value locations
    a(s2:s2:end) = a_l .* a_r + b_l .* d_r;    
    b(s2:s2:end) = a_l .* b_r + b_l .* e_r;
    c(s2:s2:end) = a_l .* c_r + b_l .* f_r + c_l;
    d(s2:s2:end) = d_l .* a_r + e_l .* d_r;
    e(s2:s2:end) = d_l .* b_r + e_l .* e_r;
    f(s2:s2:end) = d_l .* c_r + e_l .* f_r + f_l;
    
end


% do second log2(n)-1 steps
for log_s = logn-2:-1:0
    s = 2^log_s;
    s2 = 2*s;

    % pick out the right-hand values. 
    a_r = a(s2+s:s2:end);  b_r = b(s2+s:s2:end);  c_r = c(s2+s:s2:end);
    d_r = d(s2+s:s2:end);  e_r = e(s2+s:s2:end);  f_r = f(s2+s:s2:end);
    
    % and the left-hand values
    a_l = a(s2:s2:end-1);  b_l = b(s2:s2:end-1);  c_l = c(s2:s2:end-1);
    d_l = d(s2:s2:end-1);  e_l = e(s2:s2:end-1);  f_l = f(s2:s2:end-1);
    
    %Note that we need to store a copy of the RHVs because we modify them
    %in place during each iteration of the loop but need to know the
    %original values from the start of the iteration.

    % store in right hand value locations
    a(s2+s:s2:end) = a_l .* a_r + b_l .* d_r;    
    b(s2+s:s2:end) = a_l .* b_r + b_l .* e_r;
    c(s2+s:s2:end) = a_l .* c_r + b_l .* f_r + c_l;
    d(s2+s:s2:end) = d_l .* a_r + e_l .* d_r;
    e(s2+s:s2:end) = d_l .* b_r + e_l .* e_r;
    f(s2+s:s2:end) = d_l .* c_r + e_l .* f_r + f_l;
    
end


% Now we effectively have n rotation matrices, we just need to apply them
% to each point..we could loop but loops are boring!
if ~exist('ixy_other', 'var')
    xy(:,1) = a .* x_ + b .* y_ + c;
    xy(:,2) = d .* x_ + e .* y_ + f;
else
    idx = ixy_other(:,1);
    x_ = ixy_other(:,2);
    y_ = ixy_other(:,3);
    xy = [a(idx) .* x_ + b(idx) .* y_ + c(idx), ...
                    d(idx) .* x_ + e(idx) .* y_ + f(idx)];
end

% trim off the padded stuff (powers of 2, remember?)
xy(true_n+1:end,:) = [];
end


%{
Some notes on debugging...amazingly it actually didn't take long to make
this work. Here are the main things I did to debug while I was writing...

See the prefix sum diagram for 16 elements:
http://en.wikipedia.org/wiki/Prefix_sum#/media/File:Prefix_sum_16.svg
Give this function a 16-long vector and add in the following lines...
outside the loops:
    >>> dummy=1:n
inside the first loop (after s2=...)
    >>> disp([dummy(s2:s2:end); dummy(s:s2:end)]);
inside the second loop (after s2=...)
    >>> disp([dummy(s2+s:s2:end); dummy(s2:s2:end-1)]);
Note the disp stuff prints out nx2 arrays, where to top row is the destination
and the bottom row is the other operand.


You can get rot_mats and loop over point * rot_mat as follows:
>>>
    rot_mats = [a d zeros(n,1) b e zeros(n,1) c f ones(n,1)];
    rot_mats = reshape(rot_mats', [3,3,n]);
    for ii=1:n
        tmp =  rot_mats(:,:,ii) * [xy(ii,:) 1]';
        xy(ii,:) = tmp(1:2);
    end


Use dth = [0 0 0 0 0 0 0 pi/4 0 0 0 0]
Then the rot_mats should be identity up to element before pi/4,
then the same constant rotation after that.

%}