% [ a, b, c ] = fitsvd( xdata, ydata, wdata )
%   Fit a second-order polynomial to data
%
%   xdata : x-axis values
%   ydata : y-axis values
%   wdata : weights
%
% Returns coefficients a, b, c, such that
% a * x^2 + b * x + c fits to (xdata,ydata)
%
function [ a, b, c ] = fitsvd( xdata, ydata, wdata )

n = size( xdata,1 );
m = 3;

A = zeros( n, m );
b = zeros( n, 1 );

for i = 1 : n,
  A( i, 1 ) = xdata( i ) * xdata( i ) * wdata( i );
  A( i, 2 ) = xdata( i )              * wdata( i );
  A( i, 3 ) = 1.0                     * wdata( i );
  b( i )    = ydata( i )              * wdata( i );
end

[ U, S, V ] = svd( A );

a = zeros( m );
for j = 1 : m,
  for i = 1 : m,
    a( j ) = a( j ) + ( ( U( :, i )' * b ) / S( i, i ) ) * V( j, i );
  end
end

b = a( 2 );
c = a( 3 );
a = a( 1 );

return

