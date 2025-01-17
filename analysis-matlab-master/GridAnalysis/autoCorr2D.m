function ret = autoCorr2D(varargin)
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
    
    % Size of inputs
    [m,n,o] = size(in.x);

    % Set to zero where there's is no data
    x = reshape(in.x,[m*n, o]);
    x(in.nodwell(:),:) = 0;
    x = reshape(x,[m,n,o]);

    % First part: Obtain FFTs of x, the sum of squares and bins visited
    Fx = fft(fft(x,2*m-1,1),2*n-1,2);
    FsumOfSquares_x = fft(fft(x.^2,2*m-1,1),2*n-1,2);
    Fn = fft(fft(~in.nodwell,2*m-1,1),2*n-1,2);
    clear x nodwell

    % Second Part: Multiply the relevant transforms and invert to obtain the
    % equivalent convolutions
    rawCorr = fftshift_(real(ifft(ifft(Fx.*conj(Fx),[],2),[],1)));
    sums_x = fftshift_(real(ifft(ifft(bsxfun(@times,Fn,conj(Fx)),[],2),[],1)));
    sumOfSquares_x = fftshift_(real(ifft(ifft(bsxfun(@times,Fn,conj(FsumOfSquares_x)),[],2),[],1)));
    N = fftshift_(real(ifft(ifft(Fn.*conj(Fn),[],2),[],1)));
    clear Fx Fn FsumOfSquares_x

    % Account for rounding errors.
    rawCorr(abs(rawCorr)<in.tol) = 0;
    sums_x(abs(sums_x)<in.tol) = 0;
    sumOfSquares_x(sumOfSquares_x<in.tol) = 0;
    N = round(N);
    N(N<=1) = NaN;

    % Compute correlation matrix
    mapStd = sqrt(max(bsxfun(@times,sumOfSquares_x,N) - sums_x.^2,0));
    mapCovar = bsxfun(@times,rawCorr,N) - sums_x.*sums_x(2*m-1:-1:1,2*n-1:-1:1,:);

    ret.autocorrelogram = mapCovar./mapStd./mapStd(2*m-1:-1:1,2*n-1:-1:1,:);
    
    if in.PLOT_ON, PlotForDebug(ret.autocorrelogram); end;

end

function y = fftshift_(x)

y = x([ceil(end/2)+1:end 1:ceil(end/2)],[ceil(end/2)+1:end 1:ceil(end/2)],:); 

end

function PlotForDebug(ac)
    if size(ac,3) > 1
        ac = ac(:,:,1);
    end;
    im = (ac +1)/2;
    cm = jet(60);
    im = round(im*60.999 + 0.5);
    bad = repmat(isnan(ac),[1 1 3]);
    im = ind2rgb(im,cm);
    im(bad) = 1;
    image(im);
    axis image
    xlabel('bins'); ylabel('bins');

end

function in = DealWithInputs(varargin)
    defaults.x = [];
    defaults.map = 'for convenience can use .map instead of .x';
    defaults.nodwell = [];
    defaults.tol = 1e-10;
    defaults.PLOT_ON = true;
    
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
    
    % extra stuff
    if numel(in.x) == 0 && isnumeric(in.map)
        in.x = in.map; 
    end
end
