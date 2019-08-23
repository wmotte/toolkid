% [ image, nii ] = procmems( filename, outputPath, factorRO, factorPE )
% 
function [ image, nii ] = procmems( filename, outputPath, factorRO, factorPE ),

% read FID
fid = readfid( filename );

% get header
header = fid.header;
fid = fid.fid;

numberOfBlocks = header.nblocks; % number of blocks
numberOfTraces = header.ntraces; % number of traces
numberOfRO = header.np / 2; % number of points per trace
numberOfSlices = procpar( filename, 'ns' ); % number of slices
numberOfEchos = procpar( filename, 'ne' ); % number of echos
numberOfPE = procpar( filename, 'nv' ); % number of phase encoding steps
slicePositions = procpar( filename, 'pss' ); % slice positions

% zero-filling
zeroFillRO = pow2( ceil( log2( numberOfRO * factorRO ) ) );
zeroFillPE = pow2( ceil( log2( numberOfPE * factorPE ) ) );

% order of slices
[ sortedSlicePositions, sliceIndex ] = sort( slicePositions );

% order FID: RO x PE x ( slice x echo )
fid = reshape( fid, numberOfRO, numberOfEchos, numberOfSlices, numberOfPE );  % -> 256    12    19   128
fid = permute( fid, [ 1 4 3 2 ] );                                            % -> 256   128    19    12

% prepare output
image = zeros( zeroFillRO, zeroFillPE, numberOfSlices, numberOfEchos );

% transform
for slice = 1 : numberOfSlices,
  for echo = 1 : numberOfEchos,
    image( :, :, slice, echo ) = fft2d( fid( :, :, sliceIndex( slice ), echo ), zeroFillRO, zeroFillPE );
  end
end

% determine output dimensions
voxelSize = [ 1 1 1 ];

voxelSize( 1 ) = procpar( filename, 'lpe' ) / size( image, 1 ) * 10;
voxelSize( 2 ) = procpar( filename, 'lro' ) / size( image, 2 ) * 10;
sliceThickness = circshift( sortedSlicePositions, [ 0 -1 ] ) - sortedSlicePositions;
voxelSize( 3 ) = sliceThickness( 1 ) * 10;

% save to nifti
nii = make_nii( abs( image ), voxelSize, zeros( size( voxelSize ) ) );
save_nii( nii, [ outputPath '.nii' ] );

return
