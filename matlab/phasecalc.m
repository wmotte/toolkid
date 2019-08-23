% output = phasecalc( input )
%   Calculate the phase of a single input trace vector
%
function output = phasecalc( input )

phase1 = zeros( length( input ) );
phase2 = zeros( length( input ) );

for x = 1 : length( input ),
    phase1( x ) = atan2( real( input( x ) ), imag( input( x ) ) );
    phase2( x ) = mod( ( phase1( x ) + 2 * pi ), 2 * pi );
end

output = phase1;

for x = 2 : length( input ),
    dif1 = phase1( x - 1 ) - phase1( x );
    dif2 = phase2( x - 1 ) - phase2( x );
    
    if ( abs( dif1 ) < abs( dif2 ) )
        output( x ) = output( x - 1 ) + dif1;
    else
        output( x ) = output( x - 1 ) + dif2;
    end
end

output = output - output( round( length( input ) / 2 ) ) + phase1( round( length( input ) / 2 ) );

return
