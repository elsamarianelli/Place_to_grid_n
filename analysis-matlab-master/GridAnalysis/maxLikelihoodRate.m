function ret = maxLikelihoodRate(varargin)
% See https://github.com/UCL/mTint/tree/master/GridAnalysis for info
%
%%    Copyright (C) <2013>  <Ali Jeewajee> <a.jeewajee@ucl.ac.uk> &
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

    if nargin == 0, ret = DealWithInputs; return, end;
    in = DealWithInputs(varargin{:});   
    
    dwell = in.dwell;
    spike = in.spike;

    % Reshape 2D variables to make them 1D
    nVars = numel(in.combineDims);
    [n,varTrueShapes] = flattenedShape(in.combineDims,size(dwell));
    dwell = reshape(dwell,n);
    spike = reshape(spike,n);

    % Initialise each variable's vector to be ones of the relevant length
    X = arrayfun(@(a) ones(a,1),n,'UniformOutput',false);
    Xshape = diag(n-1) + 1; %1s off diagonal, n on diagonal

    % Generate perumtations of dwell for use as matricies in estimate subfunction
    dwellPerms = cell(nVars);
    for ii=1:nVars
        perm = [ii 1:nVars];
        perm(1 + ii) = [];
        dwellPerms{ii} = reshape(permute(dwell,perm),n(ii),[]);
    end

    % Generate spike subtotals for use in estimate subfunction
    spikeSubtotals = cell(nVars);
    for ii=1:nVars
        dims = 1:nVars;
        dims(ii) = [];
        spikeSubtotals{ii} = reshape(sum_(spike,dims),n(ii),1);
    end

    % Precalculate stuff for use in loop
    sumGammalnSpikePlus1 = sum_(gammaln(spike + 1),1:nVars);
    logMin = log(in.tol);
    logDwell = log(dwell);

    % Apply algorithm
    ret.converged = false;
    fit = in.fit;
    for iter = 1:in.maxIterations
        prev_fit = fit;

        % Estimate contributions of the data to different spatial correlates
        X = estimate(X, Xshape, spikeSubtotals, dwellPerms, in.tol);

        % Calculate loglikelihood of data given estimates
        expec = expected(X,Xshape,dwell);
        logExpec = logExpected(X,Xshape,logDwell,logMin);
        fit = sum_(spike.*logExpec,1:nVars) - sum_(expec,1:nVars) - sumGammalnSpikePlus1;

        % Stop if likelihood is maximised to within tolerance, otherwise cont.
        if abs(prev_fit - fit) < -in.accuracy*fit 
            fprintf(1, '\n converged, loglikelihood: %f\n', fit);
            ret.converged = true;
            break;
        end
        fprintf(1, ' %f', fit);
    end

    % If convergence is achieved rescale firing "rates" to reflect the real
    % number of spikes rather than the partial contribution to firing of a
    % given variable.
    ret.ml_rate = cell(1,nVars);
    if ~ret.converged
        fprintf(1,'\n Did not converge, loglikelihood: %f\n', fit);
    else
        total_spikes = sum_(spike,1:nVars); 

        for ii=1:nVars
            dims = 1:nVars;
            dims(ii) = [];
            pred_spikes = X{ii}' * reshape(sum_(dwell,dims),[],1); %row vector times column vector is scalar
            ret.ml_rate{ii} = X{ii}.*(total_spikes/pred_spikes);
            if numel(varTrueShapes{ii}) > 1 %if var was more than 1-d we need to restore its true shape
                ret.ml_rate{ii} = reshape(ret.ml_rate{ii},varTrueShapes{ii}); 
            end
        end

    end


end

% ######################################################################
% Estimate data given initial values.
function X = estimate(X, Xshape, spikeSubtotals, dwellPerms, tol)
    nVars = numel(X);

    for ii=1:nVars
        jj_list = 1:nVars;
        jj_list(ii) = [];
        tmp = X{jj_list(1)};
        for jj=jj_list(2:end)
            tmp = bsxfun(@times,tmp,reshape(X{jj},Xshape(jj,:)));
        end
        denom = dwellPerms{ii} * tmp(:); % matrix times vector is vector
        denom(denom<tol) = 1/tol;
        X{ii} = spikeSubtotals{ii}./denom;
    end

end

% ######################################################################
% Calculate expected values
function ret = expected(X,Xshape,dwell)
    nVars = numel(X);
    ret = 1;
    for ii = 1:nVars
        ret = bsxfun(@times,ret,reshape(X{ii},Xshape(ii,:)));
    end
    ret = ret.*dwell;
end

% ######################################################################
% Calculate log of expected values
function ret = logExpected(X,Xshape,logDwell,logMin)
    % could just do log(expected), but taking logs is slow, faster way is to
    % recalculate from scratch with sum-of-logs not log-of-products
    nVars = numel(X);
    ret = 0;
    for ii = 1:nVars
        ret = bsxfun(@plus,ret,reshape(log(X{ii}),Xshape(ii,:)));
    end
    ret = ret + logDwell;
    ret = max(ret,logMin);
end

% ######################################################################
% Sum along multiple dimensions of an nd-array
function out = sum_(in,dims)
    out = in;
    for ii=1:numel(dims)
        out = sum(out,dims(ii));
    end
end

% ######################################################################
% Calculates the shape of an nd-array when some dimensions are combined.
function [newSize,origialSizes] = flattenedShape(combineDims,oldSize)
    % combineDims is a vector of the form [1 1 2 1 ], where each element
    % specifies how many of the old dimensions should be combined to make the new.
    % Returns newSize, which can be used with reshape to apply the new shape;
    % and originalSizes, which is a cell array giving the original shapes of 
    % each of the flattened dimensions. 
    nNewDims = numel(combineDims);
    newSize = zeros(1,nNewDims);
    origialSizes = cell(1,nNewDims);

    p = 1;
    for ii=1:nNewDims
        origialSizes{ii} = oldSize(p-1 + (1:combineDims(ii)));
        newSize(ii) = prod(origialSizes{ii});
        p = p + combineDims(ii);
    end

end



function in = DealWithInputs(varargin)
    defaults.spike = [];
    defaults.dwell = [];
    defaults.maxIterations = 30;
    defaults.accuracy = 0.001;
    defaults.tol = 0.1;
    defaults.fit = 1;
    defaults.combineDims = [2 1];

    VERSION = 1;

    % Boiler plate for DealWithInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0
        in = defaults;
        return;
    end

    if isstruct(varargin{1})
        if nargin > 1 && VERSION ~= varargin{2}
            error(['%s called with version number %g, but code is version %g.\n' ...
                'Check the GitHub wiki to see what has changed, or take a risk and call ' ...
                'the function without a version number.'],mfilename,varargin{2},VERSION);
        end
        in = ModifyExistingFields(defaults,varargin{1});
    else
        in = ModifyExistingFields(defaults,varargin{:});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
