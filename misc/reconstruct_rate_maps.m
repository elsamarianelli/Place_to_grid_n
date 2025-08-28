function rate_maps = reconstruct_rate_maps(Z, pos, dims)
    K = size(Z, 2);
    rate_maps = zeros([dims, K]);
    x_edges = linspace(1, dims(2), dims(2)+1);
    y_edges = linspace(1, dims(1), dims(1)+1);

    for k = 1:K
        z_k = Z(:,k);
        valid = all(~isnan(pos),2) & ~isnan(z_k);
        row = discretize(pos(valid,2), y_edges);
        col = discretize(pos(valid,1), x_edges);
        rate_maps(:,:,k) = accumarray([row, col], z_k(valid), dims, @mean, NaN);
    end
end
