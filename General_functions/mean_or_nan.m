% Helper: mean with fallback to NaN
function val = mean_or_nan(data)
    if isempty(data)
        val = NaN;
    else
        val = mean(data);
    end
end
