function w=fixed_point_directed( u, n_max )
% fixed_point_directed(u,N) gives aproximate size distribution of weakly 
% connected components in the directed configuration model.
% Source:  "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
% CC BY 2017.

u = u / sum( u(:) );
U = u( end:-1:1, end:-1:1 );
U = U / sum( U(:) );

Ui = mpolyder( U,1 );
Uo = mpolyder( U,2 );

Ui = Ui ./ sum( Ui(:) );
Uo = Uo ./ sum( Uo(:) );

F_tol = 1e-5;

xi    = exp( -2 * pi * 1i * ( 0:(n_max) ) / (n_max+1) );

Wi = ones( 1, length( xi ) );
Wi = Wi / sum( Wi );


Wo = ones( 1, length( xi ) );
Wo = Wo / sum( Wo );
k  = 0;
nrm = 1;
while nrm > F_tol;
    
    Wi_new = xi .* mpolyval( Ui, Wo, Wi );
    Wo_new = xi .* mpolyval( Uo, Wo, Wi );
    
    k = k + 1;
    
    nrm = norm( Wi - Wi_new ) + norm( Wo - Wo_new );
    
    if mod( k, 10 )==0
        disp( nrm )
    end;
    
    
    Wo = Wo_new;
    Wi = Wi_new;

end;
%%
W = xi .* mpolyval( U, Wo, Wi ); 
w  = real( ifft( W  ) );
w(1) = [];
w(1) = 0;


