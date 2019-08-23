 
function [ D, alpha ] = DFA_main( x )

% x should be a time series of length(x) greater than 2000,and of column vector.
%A is the alpha in the paper
%D is the dimension of the time series
%n can be changed to your interest

n = 100:100:1000;

N1 = length( n );

F_n = zeros( N1, 1 );

for i=1:N1
    F_n ( i ) = DFA( x, n(i), 1 );
end

n=n';
plot( log( n ), log( F_n ) );

xlabel('n')
ylabel('F(n)')

A = polyfit( log( n( 1:end ) ), log( F_n( 1:end ) ), 1 );

alpha = A( 1 );

D = 3 - A( 1 );

return
