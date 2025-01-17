function [ f, freqs1, freqs2 ] = sft2_low( x, maxFreq, totalLen, duplicate)
% sft2_low( x, maxFreq, totalLen, duplicate)
%
% INPUTS:
% x - real valued matrix on which to perform the Fourier Transform
% maxFreq-  scalar or 2-element vector giving the max spatial freq of interest 
% totalLen- scalar or 2-element vector equivalent to the zero-padding options in
%           the standard FFT, i.e. a large number specifying the size of an 
%           equivalent full zero-padded FFT of which you only wish to compute 
%           the central section of.
% duplicate-given the rotational-complex-conjugate symmetry of the transform, in
%           most useage cases there is no need to duplicate the second half of 
%           the result.  This input is a flag controlling whether or not you 
%           require the duplication to be done.
% Note that by using 2-element vectors for maxFreq and totalLen you can specify
% different values for the two spatial dimensions.
% Note also that if maxFreq is equal to floor(size(x)/2) and x is even, you'll 
% get the nyquist frequency twice unlike in the FFT where you only get it once.
% (But in such cases you would be probably be better off using the FFT itself.)
% Note that it is possible to write a small amount of additional code to deal
% with complex valued x.
%
% RETURNS:
% f- the resulting transform (note that as with the fft you may want to 
%    fftshift this matrix though obviouslly you cannot use fftshift itself if 
%    the duplicate flag was false.)
% freqs1- a vector stating what freqs each row (i.e. dim 1) of f corresponds to 
% freqs2- a vector stating what freqs each col (i.e. dim 2) of f corresponds to 
%
% Matlab File Exchange page:
% mathworks.com/matlabcentral/fileexchange/43311-low-freq-2d-fourier-transform
%
% DM, 2013
%

[N_1,N_2] = size(x);

if numel(totalLen) == 1, totalLen(2) = totalLen(1);    end
if numel(maxFreq) == 1, maxFreq(2) = maxFreq(1);   end

k_list_1 = 0:maxFreq(1)*totalLen(1)/N_1;
k_list_2 = 0:maxFreq(2)*totalLen(2)/N_2;

f = exp(-2*pi*1i/totalLen(2) * k_list_2(:) * (0:N_2-1)) * x.'; % dim=2
f = [f ; conj(f(end:-1:2,:))];

f = f * exp(-2*pi*1i/totalLen(1) * (0:N_1-1)' * k_list_1(:)'); % dim=1
f = f.';

if duplicate
    f = [f ; conj(f(end:-1:2,[1 end:-1:2]))];
end

if nargout>1
    freqs2 = [k_list_2 -k_list_2(end:-1:2)] * N_2/totalLen(2);
    if duplicate        
        freqs1 = [k_list_1 -k_list_1(end:-1:2)] * N_1/totalLen(1);
    else
        freqs1 = [k_list_1] * N_1/totalLen(1);
    end
end

end