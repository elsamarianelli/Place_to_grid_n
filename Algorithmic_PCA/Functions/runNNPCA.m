function PC_NN = runNNPCA(NeuronxTimeMat, NumberOfPC, zero_mean)
    % Ensure reproducibility
    r_saved = NeuronxTimeMat;
    L = size(r_saved, 2);  % Number of time steps
    PC_NN = zeros(size(r_saved, 1), NumberOfPC);

    for z = 1:NumberOfPC
        disp(z)

        % Zero-mean input data
        if strcmp(zero_mean, 'spatial')
            meanInputMat = mean(r_saved, 2);  % mean over time
            InputMatFixed = r_saved - repmat(meanInputMat, 1, L);
        elseif strcmp(zero_mean, 'temporal')
            meanInputMat = mean(r_saved, 1);  % mean over neurons
            InputMatFixed = r_saved - repmat(meanInputMat, size(r_saved, 1), 1);
        else
            error('Invalid zero_mean mode. Use "spatial" or "temporal".');
        end

        % Compute non-negative principal component
        PC_NN(:, z) = NNPCA2014(InputMatFixed);

        % Break if NaNs detected
        if any(isnan(PC_NN(:, z)))
            fprintf('NaNs in PC_NN(:,%d), breaking loop...\n', z);
            break;
        end

        % Remove projection of current PC
        r_saved = InputMatFixed - PC_NN(:, z) * (PC_NN(:, z)' * InputMatFixed);
    end
end
