% y = mkpoly( a, b, c, x )
%   Evaluate second-order polynomial with coefficients a, b, c, to x-values
%
function y = mkpoly( a, b, c, x ),

y = zeros( size( x ) );
for i = 1 : size( x, 1 ),
  y( i ) = a * x( i ) * x( i ) + b * x( i ) + c;
end

return

