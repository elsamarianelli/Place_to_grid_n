function [PlaceCellsUni, PlaceCellsTanni, env,xy_field_u, xy_field_t ] = generate_place_cells(env, n_cells, dim_x, dim_y, boundary_effect, distribution_type)
% generate place cells with option for arrayed or random place cell
% distribution, returning uniformly distributed and boundary effected place
% cells. ( returned PlaceCellsTanni in case of 'array' setting just returns
% normal non warped array, need to add this in later if decide to use)

    % generate place field centres
    if strcmp(distribution_type, 'random')
        [xy_field_u, env] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, 1000);
        [xy_field_t, env] = getPlaceFieldCentres(env, n_cells, dim_x, dim_y, boundary_effect);
    elseif strcmp(distribution_type, 'array')
        [xy_field_u, env] = get_warped_array_coords(env, dim_x, dim_y, n_cells, 'uniform');
        [xy_field_t, env] = get_warped_array_coords(env, dim_x, dim_y, n_cells, 'warped');
    end

    % plot place field centres 
    figure; 
    subplot(1, 2, 1); plot(xy_field_u(:,1), xy_field_u(:,2), '.', 'MarkerSize', 15); hold on;
    subplot(1, 2, 2); plot(xy_field_t(:,1), xy_field_t(:,2), '.', 'MarkerSize', 15);

    % generate place field firing rate maps
    av_bound_dist = nanmean(env.dwmap, 'all');
    fw = fieldWidth(av_bound_dist) / 3; % field width for uniform
    PlaceCellsUni = generateUniformPCs(env, n_cells, xy_field_u, fw);
    PlaceCellsTanni = generateSanderPCs(env, n_cells, xy_field_u); 

end

