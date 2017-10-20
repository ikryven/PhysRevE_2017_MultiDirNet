function w = component_num_degenerate( u, n_max )
%component_num_degenerate( u, n_max )  aproximates the size distribution of weakly connected components in
%   directed configuration network with a degenerate degree distribution: u(k,l)=0, k>0 or u(k,l)=0, l>0.
%   u - degree distibution as a two-dimensional array, the size
%   distribution is calculated at 1:n_max.
%   Source: "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
%   Licensed under CC BY, 2017.

warning off

if ~exist('n_max')
    n_max = 1e3;
end;

F_tol = 1e-8;

%% Generating funcitons
if diff(size(u))<0
    u=u'
end;

u1 = u(2,:);
u0 = u(1,:);


U0 = u0( end:-1:1 );
U0 = U0 / sum( U0 );
U0 = U0 / sum( U0 );
U1 = u1( end:-1:1 );
U1 = U1 / sum( U1 );

xi = exp( -2 * pi * 1i * ( 0:(n_max-1) ) / (n_max) );


W1 = ones( 1, length( xi ) );
W1 = W1 / W1( end );

%% Fixed point iteration
nrm = 1;
while nrm > F_tol;
    
    W10 = W1;
    W1  = xi .* polyval( U1, W1 );
    nrm = norm( W1 - W10 );

end;

W = xi .* polyval( U0, W1 );
w  = real( ifft( W  ) );

%%
mn   = length( u ) - 1;
m0   = sum( u( 1, : ) );
m1   = sum( u( 2, : ) );
mu11 = sum( ( 0:mn )    .* u( 2, : ) ) / m1; 
mu01 = sum( ( 0:mn )    .* u( 1, : ) ) / m0; 

C = (1 + mu01/(1 - mu11) );

nn = 0:length(w)-1;
w =  w .* nn / C  ;
w( 1 ) = [];
w( 1 ) = 0;



