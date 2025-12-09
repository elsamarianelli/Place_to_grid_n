function kernel = mexican_hat_kernel(sz, sigma_exc, sigma_inh, strength)
    % sz: kernel size (e.g., 15)
    % sigma_exc: excitation Gaussian std (e.g., 1.5)
    % sigma_inh: inhibition Gaussian std (e.g., 4.0)
    % strength: scaling factor for inhibitory component (e.g., 0.8)

    [x, y] = meshgrid(linspace(-sz/2, sz/2, sz), linspace(-sz/2, sz/2, sz));
    exc = exp(-(x.^2 + y.^2) / (2 * sigma_exc^2));
    inh = exp(-(x.^2 + y.^2) / (2 * sigma_inh^2));
    kernel = exc - strength * inh;
    kernel = kernel - mean(kernel(:));  % zero-mean
end
