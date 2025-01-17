% calculate the firing rate at a given position (x, y) based on a Gaussian distribution centered at (mean_x, mean_y) with standard
% deviations sig_x and sig_y
function fr = firingRate(x, y, mean_x, mean_y, sig_x, sig_y, C)
    fr = (C / (2 * pi * sig_x * sig_y)) * exp(-(x-mean_x)^2 / (2*sig_x^2)) * exp(-(y-mean_y)^2 / (2*sig_y^2));
    % fr = exp(-(x-mean_x)^2 / (2*sig_x^2)) * exp(-(y-mean_y)^2 / (2*sig_y^2)) / (2*pi*sig_x*sig_y);
end
