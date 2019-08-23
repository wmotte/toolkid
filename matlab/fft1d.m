% transform = fft1d( fid, N )
%   Fourier transform (ifft) in dimension 1
%
%   fid : input spectrum (1D, 2D, or 3D)
%   N   : length of the spectrum in dimension 1
%         N > size( fid, 1 ): zero-filling the spectrum
%         N < size( fid, 1 ): center of the spectrum is retained
%
function transform = fft1d( fid, N );

inputSize = size( fid );
outputSize = size( fid );

% output size in the first dimension
outputSize( 1 ) = N;

% prepare output
transform = zeros( outputSize );

if N < inputSize( 1 ),
  % keep center of k-space
  length = inputSize( 1 ) - N;
  startIndex = round( length / 2 );
  endIndex = startIndex + N - 1;

  transform( :, :, : ) = fid( startIndex : endIndex, :, : );
else
  % zero-filling / copying
  transform( 1 : inputSize( 1 ), :, : ) = fid;
end

% shift, transform, center
transform = circshift( ifft( circshift( transform, round( -1 * inputSize( 1 ) / 2 ) ) ), N / 2 );

return

