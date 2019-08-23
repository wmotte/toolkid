function mask = mkmask( transform, limit ),

mask = abs( transform );
[ sortValue, sortIndex ] = sort( mask );
highIndex = floor( limit * length( sortIndex ) );
lowIndex = ceil( ( 1 - limit ) * length( sortIndex ) );
highValue = sortValue( highIndex );
lowValue = sortValue( lowIndex );

mask( mask < lowValue ) = 0;
mask( mask > highValue ) = 0;
mask( mask > 0 ) = 1;

return
