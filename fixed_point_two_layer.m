function w = fixed_point_two_layer( u, n_max )
% fixed_point_two_layer(u,N) gives aproximate size distribution of multilayer 
% connected components in the multiplex configuration model with two layers.
% Source:  "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
% CC BY 2017.

F_tol = 1e-5;

U = u( end:-1:1, end:-1:1 );
U = U / sum( U(:) );

U1 = mpolyder( U, 1 );
U2 = mpolyder( U, 2 );

U1 = U1 ./ sum( U1(:) );
U2 = U2 ./ sum( U2(:) );

xi = exp( -2 * pi * 1i * ( 0:(n_max) ) / (n_max+1) );

W1 = ones( 1, length( xi ) );
W1 = W1 / sum( W1 );

W2 = ones( 1, length( xi ) );
W2 = W2 / sum( W2 );

k = 0;
nrm = 1;
while nrm > F_tol;
    
    W1_new = xi .* mpolyval( U1, W1, W2 );
    W2_new = xi .* mpolyval( U2, W1, W2 );
    
    k = k + 1;
    
    nrm = norm( W1 - W1_new ) + norm( W2 - W2_new );
    
    if mod( k, 10 ) == 0
        disp( nrm )
    end;
    
    W1 = W1_new;
    W2 = W2_new;

end;
%%

U = u( end : - 1 : 1, end : - 1 : 1 );
W = xi .* mpolyval( U, W1, W2 ); 
w  = real( ifft( W  ) );
w  =  w( 2:end );
w(1) = 0;
