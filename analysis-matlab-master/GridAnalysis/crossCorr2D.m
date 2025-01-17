function ret = crossCorr2D(varargin)
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
    [mx,nx,ox] = size(in.x);
    [my,ny,oy] = size(in.y);

    % Set to zero where there's is no data
    x = reshape(in.x,[mx*nx, ox]);
    x(in.nodwell_x(:),:) = 0;
    x = reshape(x,[mx,nx,ox]);

    y = reshape(in.y,[my*ny, oy]);
    y(in.nodwell_y(:),:) = 0;
    y = reshape(y,[my,ny,oy]);

    % First part: Obtain FFTs of x, the sum of squares and bins visited
    Fx = fft(fft(x,2*my-1,1),2*ny-1,2);
    FsumOfSquares_x = fft(fft(x.^2,2*my-1,1),2*ny-1,2);
    Fnx = fft(fft(~in.nodwell_x,2*my-1,1),2*ny-1,2);
    clear x nodwell_x

    Fy = fft(fft(y,2*mx-1,1),2*nx-1,2);
    FsumOfSquares_y = fft(fft(y.^2,2*mx-1,1),2*nx-1,2);
    Fny = fft(fft(~in.nodwell_y,2*mx-1,1),2*nx-1,2);
    clear y nodwell_y

    % Second Part: Multiply the relevant transforms and invert to obtain the
    % equivalent convolutions
    rawCorr = fftshift_(real(ifft(ifft(Fx.*conj(Fy),[],2),[],1)));
    sums_x = fftshift_(real(ifft(ifft(bsxfun(@times,Fx,conj(Fny)),[],2),[],1)));
    sums_y = fftshift_(real(ifft(ifft(bsxfun(@times,Fnx,conj(Fy)),[],2),[],1)));
    sumOfSquares_x = fftshift_(real(ifft(ifft(bsxfun(@times,FsumOfSquares_x,conj(Fny)),[],2),[],1)));
    sumOfSquares_y = fftshift_(real(ifft(ifft(bsxfun(@times,Fnx,conj(FsumOfSquares_y)),[],2),[],1)));
    N = fftshift_(real(ifft(ifft(Fnx.*conj(Fny),[],2),[],1)));
    clear Fx Fy Fnx Fny FsumOfSquares_x FsumOfSquares_y

    % Account for rounding errors.
    rawCorr(abs(rawCorr)<in.tol) = 0;
    sums_x(abs(sums_x)<in.tol) = 0;
    sums_y(abs(sums_y)<in.tol) = 0;
    sumOfSquares_x(abs(sumOfSquares_x)<in.tol) = 0;
    sumOfSquares_y(abs(sumOfSquares_y)<in.tol) = 0;
    N = round(N);
    N(N<=1) = NaN;

    % Compute correlation matrix
    mapStd_x = sqrt(bsxfun(@times,sumOfSquares_x,N) - sums_x.^2);
    mapStd_y = sqrt(bsxfun(@times,sumOfSquares_y,N) - sums_y.^2);
    mapCovar = bsxfun(@times,rawCorr,N) - sums_x.*sums_y;

    ret.crosscorrelogram = mapCovar./(mapStd_x.*mapStd_y);

    if in.PLOT_ON, PlotForDebug(ret.crosscorrelogram); end;

end

function y = fftshift_(x)
    y = x([ceil(end/2)+1:end 1:ceil(end/2)],[ceil(end/2)+1:end 1:ceil(end/2)],:); 
end


function PlotForDebug(cc)
    if size(cc,3) > 1
        cc = cc(:,:,1);
    end;
    im = (cc +1)/2;
    cm = jet(60);
    im = round(im*60.999 + 0.5);
    bad = repmat(isnan(cc),[1 1 3]);
    im = ind2rgb(im,cm);
    im(bad) = 1;
    image(im);
    axis image
    xlabel('bins'); ylabel('bins');
end

function in = DealWithInputs(varargin)
    defaults.x = [];
    defaults.y = [];
    defaults.nodwell_x = [];
    defaults.nodwell_y = [];
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

end

