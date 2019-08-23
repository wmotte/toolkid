% fid = readfid( filename )
%   load Varian FID file
%
%   filename : path to FID folder (without .fid)
%   output : struct( 'header', header, 'fid', fid, 'dim', dim )
%
function fid = readfid( filename )

% open file
id = fopen( [ filename '.fid/fid' ], 'r', 'b' );

% read FID header
header = struct();
header.nblocks   = fread( id, 1, 'int32' ); % number of blocks in file
header.ntraces   = fread( id, 1, 'int32' ); % number of traces per block
header.np        = fread( id, 1, 'int32' ); % number of elements per trace
header.ebytes    = fread( id, 1, 'int32' ); % number of bytes per element
header.tbytes    = fread( id, 1, 'int32' ); % number of bytes per trace
header.bbytes    = fread( id, 1, 'int32' ); % number of bytes per block
header.vers_id   = fread( id, 1, 'int16' ); % software version and file_id status bits
header.status    = fread( id, 1, 'int16' ); % status of whole file
header.nbheaders = fread( id, 1, 'int32' ); % number of block headers

% matrix of np/2 complex points in each trace and block
dim = [ header.np / 2, header.ntraces, header.nblocks ];

% initialize FID
fid = zeros( dim );

% either 16-bit or 32-bit acquisition
precision = 'float32';
if header.ebytes == 2,
  precision = 'int16';
end

% length of block
blocksize = header.ntraces * header.np;

% read all blocks
numberOfPoints = header.np / 2;
for i = 1 : header.nblocks,
  blockhead = fread( id, 7, 'int32' ); % skip header

  for j = 1 : header.ntraces,
    trace = reshape( fread( id, header.np, precision ), 2, numberOfPoints );
    fid( :, j, i ) = complex( trace( 1, : ), trace( 2, : ) );
  end
end

% close
fclose( id );

% return as struct
fid = struct( 'header', header, 'dim', dim, 'fid', fid );

return

