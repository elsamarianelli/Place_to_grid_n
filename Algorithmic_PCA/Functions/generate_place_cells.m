function [PlaceCells, In] = generate_place_cells(In, pc_density, pc_size)
% GENERATE_PLACE_CELLS Generates place cells with controllable density and field size
% 
%   [PlaceCells, env] = generate_place_cells(In, pc_density, pc_size)
%
%   This function generates place cells in an environment, allowing:
%     - Control over whether the place cell distribution varies with distance to boundaries (pc_density)
%     - Control over whether the place field size varies with boundary distance (pc_size)
%
%   INPUTS:
%     In          - structure with environment and cell parameters
%     pc_density  - boolean, true if place field density should vary near boundaries
%     pc_size     - boolean, true if place field size should vary near boundaries
%
%   OUTPUTS:
%     PlaceCells  - generated place cell structure
%     In.env      - updated environment structure

    % [1] Set place field centres 
    if pc_density
        % Boundary-modulated place field density
        [xy_field, In.env] = getPlaceFieldCentres(In.env, In.n_cells, In.dim_x, In.dim_y, In.bound_ctrl);
    else
        % Uniform place field density
        [xy_field, In.env] = getPlaceFieldCentres(In.env, In.n_cells, In.dim_x, In.dim_y, 1000);
    end

    % [2] Set place field widths 
    if pc_size
        % Boundary-modulated place field size
        PlaceCells = generateSanderPCs(In.env, In.n_cells, xy_field);
    else
        % Uniform place field size
        av_bound_dist = nanmean(In.env.dwmap, 'all');
        fw = fieldWidth(av_bound_dist) / In.pf_width_cntrl;
        PlaceCells = generateUniformPCs(In.env, In.n_cells, xy_field, fw);
    end

    %  % [Optional] Plot field centres 
    % figure;
    % scatter(xy_field(:,1), xy_field(:,2), 20, 'filled');
    % title('Place Field Centres');
    % axis equal tight;

end

