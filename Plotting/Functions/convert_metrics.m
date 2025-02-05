function [scale, orientation, gridness, ellipticity] = convert_metrics(metrics)
% changing format of metrics to make easier for plotting, can take more
% metrics in output if needed

num_boxes = size(metrics, 2);

 
example_idx = find(cellfun(@isstruct, metrics), 1);
if ~isempty(example_idx)
    sample_fields = fieldnames(metrics{example_idx}); 
else
    error('No valid structures found in metrics.');
end

metrics_extracted = struct();

for field_idx = 1:numel(sample_fields)
    metric_name = sample_fields{field_idx}; 
    metrics_extracted.(metric_name) = cell(1, num_boxes);

    % Loop through each box (1 to 9)
    for box_idx = 1:num_boxes
        valid_rows = find(cellfun(@isstruct, metrics(:, box_idx)));
        
        metric_values = arrayfun(@(r) metrics{r, box_idx}.(metric_name), valid_rows, 'UniformOutput', false);
        
        metrics_extracted.(metric_name){box_idx} = cell2mat(metric_values); 
    end
end

scale = metrics_extracted.scale;
orientation = metrics_extracted.orientation;
gridness = metrics_extracted.gridness;
ellipticity = metrics_extracted.ellipticity;

end