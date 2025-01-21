function fmap_avg = get_fmap_avg(PlaceCells)
% get mean firing rate map across all place cells

[m, n] = size(PlaceCells{1, 1}.fmap);
fmap_sum = zeros(m, n);
for i = 1:200
    fmap_sum = fmap_sum + PlaceCells{1, i}.fmap;
end
fmap_avg = fmap_sum / 200;
figure;
imagesc(fmap_avg)

end