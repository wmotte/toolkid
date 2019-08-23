% transform = phasewithmap( trace, phasemap )
%   Apply a phasemap to a trace vector
%
function transform = phasewithmap( trace, phasemap );

transform = ( cos( phasemap ) .* real( trace ) + sin( phasemap ) .* imag( trace ) ) + i *( cos( phasemap ) .* imag( trace ) - sin( phasemap ) .* real( trace ) );

