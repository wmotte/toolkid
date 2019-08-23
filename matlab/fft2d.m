% transform = fft2d( fid, N, M )
%   2-D Fourier transform
%
function transform = fft2d( fid, N, M ),

transform = fft1d( fft1d( fid, N )', M )';

return

