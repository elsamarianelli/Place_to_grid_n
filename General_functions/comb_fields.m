function [map] = comb_fields(maps,combo)

map = zeros(size(maps,[1 2]));

for ii = 1:size(maps,3)
    map = map + combo(ii)*maps(:,:,ii);
end

map = map./ii;
