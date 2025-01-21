function [ellipticity, fit_result] = assess_ellipticity(closestPeaksCoord, ax, plotting)
% Extract coordinates
x = closestPeaksCoord(:, 1);
y = closestPeaksCoord(:, 2);

% Fit an ellipse using Direct Least Squares Method
% Add 'fit_ellipse' function available from MATLAB File Exchange or custom implementation

if strcmp(plotting, 'yes')
    fit_result = fit_ellipse(x, y, ax); % Fit ellipse to the data points
else
    fit_result = fit_ellipse(x, y); % Fit ellipse to the data points
end
% Calculate ellipticity
ellipticity = fit_result.a / fit_result.b; % Ratio of semi-major to semi-minor axis

% hold on
% plot(ax, x, y, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Original points

% % Display results
% fprintf('Semi-Major Axis (a): %.2f\n', a);
% fprintf('Semi-Minor Axis (b): %.2f\n', b);
% fprintf('Ellipticity: %.2f\n', ellipticity);
% fprintf('Orientation (phi, degrees): %.2f\n', rad2deg(phi));
end