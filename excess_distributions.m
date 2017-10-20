
[ Y X ] = meshgrid( 0:( size( u, 2 ) -1 ), 0:( size( u, 1 ) - 1 ) );
u( 1, 1 ) = 0;
u = u / sum( u(:) );

[MU   SU ] = get_S_mu( u );

[ ma10   Sa10 ] = get_S_mu( X .* u / MU( 1 ) );
[ ma01   Sa01 ] = get_S_mu( Y .* u / MU( 1 ) );

m10 = sum( sum( X          .*u ) ); 
m20 = sum( sum( X.^2       .*u ) ); 
m02 = sum( sum( Y.^2       .*u ) ); 
m11 = sum( sum( X .* Y     .*u ) ); 
m21 = sum( sum( X.^2 .*Y   .*u ) );
m12 = sum( sum( X .* Y.^2 .* u ) );

