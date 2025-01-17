% generate environment
function env = generate_environment(dim_x, dim_y)
    n_polys = 1;
    polys = cell(n_polys, 1);
    polys{1} = [0 0, (dim_x-2) 0, (dim_x-2) (dim_y-2), 0 (dim_y-2), 0 0] + 2;
    env = GenerateEnv(polys, dim_x, dim_y);
end
