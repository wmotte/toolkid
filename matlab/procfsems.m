% [ image, nii ] = procfsems( filename, outputPath, factorRO, factorPE )
%
function [ image, nii ] = procfsems( filename, outputPath, factorRO, factorPE ),

% read FID
fid = readfid( filename );

% get header
header = fid.header;
fid = fid.fid;

numberOfBlocks = header.nblocks; % number of blocks
numberOfTraces = header.ntraces; % number of traces
numberOfRO = header.np / 2; % number of points per trace
numberOfSlices = procpar( filename, 'ns' ); % number of slices
numberOfEchos = procpar( filename, 'etl' ); % echo train length
numberOfPE = numberOfTraces / numberOfSlices; % number of phase encoding steps
numberOfEchoTraces = numberOfPE / numberOfEchos; % number of phase encoding steps per echo
slicePositions = procpar( filename, 'pss' ); % slice positions

% phase-encoding table
phaseTable = procpar( filename, 'pelist' );
phaseTable = reshape( phaseTable, numberOfEchos, numberOfEchoTraces );
minimumPhaseTableIndex = min( min( phaseTable ) ) - 1;

% zero-filling
zeroFillRO = pow2( ceil( log2( numberOfRO * factorRO ) ) );
zeroFillPE = pow2( ceil( log2( numberOfPE * factorPE ) ) );

% order of slices
[ sortedSlicePositions, sliceIndex ] = sort( slicePositions );

% order FID: RO x ( echo x slices x phase ) x blocks
fid = reshape( fid, numberOfRO, numberOfEchos, numberOfSlices, numberOfEchoTraces, numberOfBlocks );

% prepare output
image = zeros( zeroFillRO, zeroFillPE, numberOfSlices, numberOfBlocks );

% transform
for block = 1 : numberOfBlocks, 
  
  % re-order phase lines
  fidBlock = zeros( numberOfRO, numberOfSlices, numberOfPE );

  for echo = 1 : numberOfEchos,
    for phase = 1 : numberOfEchoTraces,
      % look-up
      row = phaseTable( echo, phase ) - minimumPhaseTableIndex;
 
      % fill block FID 
      fidBlock( :, :, row ) = squeeze( fid( :, echo, :, phase, : ) );
    end
  end

  % re-order: RO x PE x slice
  fidBlock = permute( fidBlock, [ 1 3 2 ] );

  % transform
  for slice = 1 : numberOfSlices,
    image( :, :, slice, block ) = fft2d( fidBlock( :, :, sliceIndex( slice ) ), zeroFillRO, zeroFillPE ); 
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
