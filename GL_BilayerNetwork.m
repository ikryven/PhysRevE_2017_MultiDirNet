function w = GL_BilayerNetwork( u, N )
% GL_BilayerNetwork( u, N ) computes exact size distribution in the multiplex
% configuration model with two layers.
% Source:  "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
% CC BY 2017.

u( N + 2, N + 2 ) = 0;

[ Y X ] = meshgrid( 0:N+1 );

Ux = X .* u;
Uy = Y .* u;

Ux = Ux / sum( Ux(:) );
Uy = Uy / sum( Uy(:) );

Ux = Ux( 2:end, :        );
Uy = Uy(        :, 2:end );

Ux = Ux( 1:N, 1:N );
Uy = Uy( 1:N, 1:N );

[Y X] = meshgrid( 0:N-1 );

FUy  = fft2( Uy );
FUx  = fft2( Ux );

FXUx = fft2( X .* Ux );
FYUy = fft2( Y .* Uy );

FYUx = fft2( Y .* Ux );
FXUy = fft2( X .* Uy );

FU   = fft2( u( 1:N, 1:N ) );

D =  ( FUx - FXUx ) .* ( FUy - FYUy ) - FXUy .* FYUx;
w = u( 1, 1 );
for n = 2:N-1;
    n
    w( n+1 ) = 0;
    for i = 0:n-1
        j = n - i - 1;
        
        tmp = ifft2( FU .* D .* FUx .^ ( i - 1 ) .* FUy .^ ( j - 1 )   );        
        w( n ) = w( n )  + real( tmp( i  +1, j  + 1 ) );
        
    end;
end


 