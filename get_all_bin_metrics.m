function [metrics_nine, metrics_three] = get_all_bin_metrics(map, shape)
    % Compute SAC once
    [sac, metrics_nine] = get_binned_metrics(map, shape, 'nine');

    % Extract metrics for 'three' from 'nine'
    corners = {[1,1], [1,3], [3,1], [3,3]};
    edges   = {[1,2], [2,1], [2,3], [3,2]};
    center  = [2,2];

    fields = fieldnames(metrics_nine{1,1});
    for f = 1:numel(fields)
        fn = fields{f};

        % Collect corner metrics
        c_vals = cellfun(@(idx) metrics_nine{idx(1), idx(2)}.(fn), corners, 'UniformOutput', false);
        c_vals = cell2mat(c_vals(~cellfun(@isempty, c_vals)));

        % Collect edge metrics
        e_vals = cellfun(@(idx) metrics_nine{idx(1), idx(2)}.(fn), edges, 'UniformOutput', false);
        e_vals = cell2mat(e_vals(~cellfun(@isempty, e_vals)));

        % Center metric
        center_val = metrics_nine{center(1), center(2)}.(fn);

        % Assign to output struct
        metrics_three(1).(fn) = mean(c_vals, 'omitnan');
        metrics_three(2).(fn) = mean(e_vals, 'omitnan');
        metrics_three(3).(fn) = center_val;
    end
end
