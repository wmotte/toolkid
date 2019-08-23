% phasemap = epiphasemap( trace, zeroFillRO )
%   Build a phase map for EPI phase correction
%
%   trace : FID trace (RO x PE matrix)
%   zeroFillRO : output vector length in RO-direction
%
function phasemap = epiphasemap( trace, zeroFillRO ),

numberOfRO = size( trace, 1 );
numberOfPE = size( trace, 2 );

% transform in RO direction
transform = fft1d( trace, zeroFillRO );

% mask background
mask = mkmask( transform, 0.95 );

% calculate fitting weights
weight = abs( transform ).^2;
weight = weight / max( weight );

% determine full x-axis
x = [ 1 : zeroFillRO ]';

% initialize phasemap
phasemap = zeros( zeroFillRO, numberOfPE );

for i = 1 : numberOfPE,
  % calculate current phase line's phase
  phase = phasecalc( transform( :, i ) );

  % fit a second-order polynomial a*x^2 + b*x + c
  maskX = find( mask( :, i ) > 0 );
  [ a, b, c ] = fitsvd( x( maskX ), phase( maskX ), weight( maskX ) );
  phasemap( :, i ) = mkpoly( a, b, c, x );

  % correct with mean phase
  meanphase = phasecalc( phasewithmap( transform( :, i ), phasemap( :, i ) ) );

  % store in phasemap
  phasemap( :, i ) = phasemap( :, i ) - meanphase( round( zeroFillRO / 2 ) );
end

return

