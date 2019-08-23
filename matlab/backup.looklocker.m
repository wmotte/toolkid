% Name : procsems_ll
% Optimized by Wim (26-05-2009)
%
% RUN as: looklocker( 'gems_01', 'processed/gems', 1, 1 )
% -------------------------------------------------------
function [im] = looklocker( filename, exportname, zffactor1, zffactor2 );

fid                 = read_fid( filename );                   % size( fid ) -> 64 38400
ne                  = read_procpar( filename, 'ne' );         % 24
nv                  = read_procpar( filename, 'nv' );         % 64      phase-encoding steps
np                  = read_procpar( filename, 'np' ) * 0.5;   % 64      readout steps: is always RO * 2 ! (so devide by 2).
ns                  = read_procpar( filename, 'ns' );         % 25      slices
pss                 = read_procpar( filename, 'pss' );        % -0.15 -0.05 0.05 0.15 0.25   [...] 1    (25 values...)
[pssord, sortindex] = sort( pss );                            % sortindex ->  1 14 2 15 3 16 [...] 13   (25 values...)

% reshape...
fid = reshape( fid, nv, ns, ne, np );         % -> 64    25    24    64
fid = permute( fid, [1,4,2,3]      );         % -> 64    64    25    24
fid = reshape( fid, nv, np, ns, ne );         % -> 64    64    25    24

% zerofill...
zf1 = pow2( ceil( log2( nv ) ) ) * zffactor1;
zf2 = pow2( ceil( log2( np ) ) ) * zffactor2; 

im  = zeros( zf1, zf2, ns, ne );

% fft...
for slicecount=1:ns;       % 1:25
    for arraycount=1:ne;   % 1:24
      im( :, :, slicecount, arraycount ) = fft2d( fid ( :, :, sortindex( slicecount ), arraycount ), zf1, zf2 );
    end;
end;

% save...
exportbfloat( abs( im ), exportname );

