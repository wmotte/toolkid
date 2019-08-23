% episortro( fid, filename )
%   Re-order traces and points in an EPI FID
%
%   fid : FID (np x ntraces x nblocks)
%   filename : path to the FID (for procpar)
%
function [ fid, np, nv ] = episortro( fid, filename ),

% get header
header = fid.header;
dim = fid.dim;
fid = fid.fid;

% load RO-pattern table
patternRO = procpar( filename, 'ropat' );
load( [ '/home1/maurits/vnmrsys/epi_indices/' patternRO '.mat' ], 'EPI_indices', 'nv', 'np' );

numberOfRO = np; % number of points per trace
numberOfPE = nv; % number of phase encoding steps
numberOfBlocks = header.nblocks; % number of blocks
numberOfTraces = header.ntraces; % number of traces

% order FID
orderFid = zeros( numberOfRO, numberOfPE, numberOfTraces, numberOfBlocks );

for i = 1 : numberOfRO,
  for j = 1 : numberOfPE,
    orderFid( i, j, :, : ) = fid( EPI_indices( i, j ), :, : );
  end
end

fid = struct( 'header', header, 'dim', dim, 'fid', orderFid );

return

