function val = get_field(metrics, fname)
% Safe retrieval of metric fields
val = NaN;
if isfield(metrics, fname) && ~isempty(metrics.(fname))
    val = metrics.(fname);
end
end