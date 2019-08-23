% [ image, nii ] = proc3d( filename, outputPath, factorRO, factorPE, factorPE2 )
%
function [ image, nii ] = proc3d( filename, outputPath, factorRO, factorPE, factorPE2 ),

% read FID
fid = readfid( filename );

% get header
header = fid.header;
fid = fid.fid;

% number of blocks
numberOfBlocks = header.nblocks;

% number of traces
numberOfTraces = header.ntraces;

% number of points per trace
numberOfRO = header.np / 2;

% number of slices
numberOfPE2 = procpar( filename, 'nv2' );

% number of phase encoding steps
numberOfPE = procpar( filename, 'nv' );

% zero-filling
zeroFillRO  = pow2( ceil( log2( numberOfRO  * factorRO ) ) );
zeroFillPE  = pow2( ceil( log2( numberOfPE  * factorPE ) ) );
zeroFillPE2 = pow2( ceil( log2( numberOfPE2 * factorPE2 ) ) );

% re-order FID: RO x PE x PE2 x blocks
fid = reshape( fid, numberOfRO, numberOfPE, numberOfPE2, numberOfBlocks );

% prepare output
image = zeros( zeroFillRO, zeroFillPE, zeroFillPE2, numberOfBlocks );

% transform
for block = 1 : numberOfBlocks, 

  image( :, :, :, block ) = fft3d( fid( :, :, :, numberOfBlocks ), zeroFillRO, zeroFillPE, zeroFillPE2 );

end

% determine output dimensions
voxelSize = [ 1 1 1 ];

voxelSize( 1 ) = procpar( inputPath, 'lpe' ) / size( image, 1 ) * 10;
voxelSize( 2 ) = procpar( inputPath, 'lro' ) / size( image, 2 ) * 10;
voxelSize( 3 ) = procpar( inputPath, 'lpe2' ) / size( image, 3 ) * 10;

% save to nifti
nii = make_nii( abs( image ), voxelSize, zeros( size( voxelSize ) ) );
save_nii( nii, [ outputPath '.nii' ] );

return
