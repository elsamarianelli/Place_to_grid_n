function rate = lookup_ideal_rate(xy, scale, orientation, phase_xy)
% See readme.md in this folder for details, or view it on github with all
% the nice formatting and images.
%
% DM, May 2015.
%

persistent unit_grid
nBins = 100; % this is the number of bins along each axis used for rseolution
             % within a single grid unit. Note that although the unit_grid is
             % nBins x nBins, it should not be interpreted as a square, but as a
             % parallelogram...this is what the M transformation acheives
             % (We use it on the xy data and use the inverse on the
             % unit_grid xy data.)

% first do inverse pahse offset
xy_ = bsxfun(@minus, xy, phase_xy*scale); 

% Then do inverse rotation, followed by projection down to unit tile
R = [cos(-orientation),  -sin(-orientation);
     sin(-orientation) cos(-orientation)];
M = [1 0; -1/sqrt(3) 2/sqrt(3)];
bin_idx = mod(round(xy_*(R*M) * (nBins / scale)), nBins) +1; % and convert to bins

if isempty(unit_grid)
    % Make unit grid with 0-orient, 
    [xx,yy] = meshgrid(linspace(0,1,nBins));
    xy_mat = [xx(:) yy(:)];
    xy_mat = xy_mat/M;
    omegaRadius = (1/(cosd(30))) * 2 * pi; %Convert to grid scale and get appropriate vector
    omegaTheta = [-pi/6, pi/6, 3*pi/6]; %Angle of bands relative to one another
    [omegaX, omegaY] = pol2cart(omegaTheta, omegaRadius); %Vector defining pitch and direction of each band
    omega = [omegaX;omegaY]'; % [3x2] so three xy pairs
    temp = cos(omega*xy_mat'); % [3 x nBins] each row being a single band
    temp = sum(temp); %Sum accross bands
    temp = temp.*(temp>0);
    unit_grid = temp(:);
    unit_grid = unit_grid / mean(unit_grid); 
end

rate = unit_grid(bin_idx(:,1) + nBins*(bin_idx(:,2) -1)); % inlined "sub2ind"

end
