% [ image, nii ] = procepi( referencePath, inputPath, outputPath, factorRO, factorPE )
%
function [ image, nii ] = procepi( referencePath, inputPath, outputPath, factorRO, factorPE ),

% read reference
referenceFid = readfid( referencePath );

% re-order data points
[ referenceFid, np, nv ] = episortro( referenceFid, referencePath );

referenceHeader = referenceFid.header;
referenceFid = referenceFid.fid;

numberOfTraces = referenceHeader.ntraces; % number of traces
numberOfSlices = procpar( referencePath, 'ns' ); % number of slices
numberOfRO = np; % number of points per trace
numberOfPE = nv; % number of phase encoding steps

% zero-filling
zeroFillRO = pow2( ceil( log2( numberOfRO * factorRO ) ) );
zeroFillPE = pow2( ceil( log2( numberOfPE * factorPE ) ) );

slicePositions = procpar( inputPath, 'pss' ); % slice positions

% order of slices
[ sortedSlicePositions, sliceIndex ] = sort( slicePositions );

% initialize phasemap
phasemap = zeros( zeroFillRO, numberOfPE, numberOfTraces );

% build phasemap
for i = 1 : numberOfTraces,
  phasemap( :, :, i ) = epiphasemap( referenceFid( :, :, sliceIndex( i ) ), zeroFillRO );
end

% read FID
fid = readfid( inputPath );

% re-order data points
[ fid, np, nv ] = episortro( fid, inputPath );

% get header
header = fid.header;
fid = fid.fid;

numberOfBlocks = header.nblocks; % number of blocks

% prepare output
image = zeros( zeroFillRO, zeroFillPE, numberOfSlices, numberOfBlocks );

% transform
for block = 1 : numberOfBlocks, 
  for slice = 1 : numberOfSlices,
    z = sliceIndex( slice );

    plane = zeros( zeroFillRO, numberOfPE );

    for phase = 1 : numberOfPE,
      trace = fft1d( squeeze( fid( :, phase, z, block ) ), zeroFillRO );
      plane( :, phase ) = phasewithmap( trace, phasemap( :, phase, slice ) );
    end

    image( :, :, slice, block ) = fft1d( plane', zeroFillPE )';
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

