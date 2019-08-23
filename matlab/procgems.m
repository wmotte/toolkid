% [ image, nii ] = procgems( filename, outputPath, factorRO, factorPE )
%
function [ image, nii ] = procgems( filename, outputPath, factorRO, factorPE ),

% read FID
fid = readfid( filename );

% get header
header = fid.header;
fid = fid.fid;

numberOfBlocks = header.nblocks; % number of blocks
numberOfTraces = header.ntraces; % number of traces
numberOfRO = header.np / 2; % number of points per trace
numberOfSlices = procpar( filename, 'ns' ); % number of slices
numberOfPE = numberOfTraces / numberOfSlices; % number of phase encoding steps
slicePositions = procpar( filename, 'pss' ); % slice positions

% zero-filling
zeroFillRO = pow2( ceil( log2( numberOfRO * factorRO ) ) );
zeroFillPE = pow2( ceil( log2( numberOfPE * factorPE ) ) );

% re-order FID: RO x PE x slices x blocks
fid = reshape( fid, numberOfRO, numberOfSlices, numberOfPE, numberOfBlocks );
fid = permute( fid, [ 1, 3, 2, 4 ] );

% order of slices
[ sortedSlicePositions, sliceIndex ] = sort( slicePositions );

% prepare output
image = zeros( zeroFillRO, zeroFillPE, numberOfSlices, numberOfBlocks );

% transform
for block = 1 : numberOfBlocks, 
  for slice = 1 : numberOfSlices,

    image( :, :, slice, block ) = fft2d( fid( :, :, sliceIndex( slice ), numberOfBlocks ), zeroFillRO, zeroFillPE );

  end
end

% determine output dimensions
voxelSize = [ 1 1 1 ];

voxelSize( 1 ) = procpar( inputPath, 'lpe' ) / size( image, 1 ) * 10;
voxelSize( 2 ) = procpar( inputPath, 'lro' ) / size( image, 2 ) * 10;
sliceThickness = circshift( sortedSlicePositions, [ 0 -1 ] ) - sortedSlicePositions;
voxelSize( 3 ) = sliceThickness( 1 ) * 10;

% save to nifti
nii = make_nii( abs( image ), voxelSize, zeros( size( voxelSize ) ) );
save_nii( nii, [ outputPath '.nii' ] );

return
