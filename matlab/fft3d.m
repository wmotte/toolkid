% transform = fft3d( fid, N, M, P )
%   3-D Fourier transform
%
function transform = fft3d( fid, N, M, P ),

inputSize = size( fid );
outputSize = [ N M P ];

% prepare output
transform = zeros( outputSize );

% for each PE2 line, 2-D transform
for i = 1 : inputSize( 3 ),
  transform( :, :, i ) = fft2d( fid( :, :, i ), N, M );
end

% permute, transform in PE2-direction, permute back
transform = permute( fft1d( permute( transform, [ 3 1 2 ] ), P ), [ 2 3 1 ] );

return

