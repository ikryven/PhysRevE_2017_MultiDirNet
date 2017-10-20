function w = components_lagrange_DD( u, n_max )
%components_lagrange_DD( u, n_max ) gives  exact size distribution of weakly connected components in
%   directed configuration networks with a degenerate degree distribution.
%   u - degree distibution as a two-dimensional array, the size
%   distribution is calculated at 1:n_max.
%   Source: "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
%   Licensed under CC BY, 2017.

    %% Initialize
    
    mn   = length( u ) - 1;
    m0   = sum( u( 1, : ) );
    m1   = sum( u( 2, : ) );

    mu11 = sum( ( 0:mn )    .* u( 2, : ) ) / m1; 
    mu01 = sum( ( 0:mn )    .* u( 1, : ) ) / m0; 
  
    C = (1 + mu01/(1 - mu11) );
    
    u0 = u(1,:)/sum(u(1,:));
    u1 = u(2,:)/sum(u(2,:));

    u0( n_max + 2 : end ) = [ ];
    u0( n_max + 1 ) = 0;

    u1( n_max + 2 : end ) = [ ];
    u1( n_max + 1 ) = 0;

    
    %% Convolution powers
    nn=0:length(u1)-1;
    fft_u1   = fft( u1 );
    fft_ku0  = fft( nn .* u0 );
    fft_u1n = 1;
    
    w  = zeros( 1, n_max );
    for n = 2 : n_max - 1; 
        
        fft_u1n = fft_u1n .* fft_u1;
        u1n     = real( ifft( fft_u1n .* fft_ku0 ) );
        w( n ) = n / ( n - 1 ) * u1n( n ) ;
           
    end;
    
    w = w/C;
      

    
    