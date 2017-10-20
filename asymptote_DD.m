function wa = asymptote_DD( u, n )
% wa = asymptote_DD( u, n ) computes the asymptote of the size distribution of weakly connected components in
% the directed configuration model featuring degenerate degree distibution: u(k,l)=0, k>0 or u(k,l)=0, l>0.
% u - degree distibution as a two-dimensional array.
% n - one-dimensional array of points where the asymptote is to be calculated.
% Source:  "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
% CC BY 2017. 

if diff(size(u))<0
    u = u'
end;

mn   = length( u ) - 1;
m0   = sum( u( 1, : ) );
m1   = sum( u( 2, : ) );

mu12 = sum( ( 0:mn ).^2 .* u( 2, : ) ) / m1; 
mu02 = sum( ( 0:mn ).^2 .* u( 1, : ) ) / m0; 

mu11 = sum( ( 0:mn )    .* u( 2, : ) ) / m1; 
mu01 = sum( ( 0:mn )    .* u( 1, : ) ) / m0; 

E1 = ( mu11 - 1 )^2 / ( 2 * ( mu12-mu11^2 ) );
L0 = ( mu01 * ( -1 + mu11 ) ) / ( ( mu11 - mu01 -1 ) * sqrt( 2 * pi*( -mu11^2 + mu12 ) ) );
E0 = ( mu01 +  mu02 - mu01 * mu11  ) * (  mu11 - 1 )  / ( mu01 * ( mu12 - mu11^2   ) );

wa  = L0 * n.^( -1/2 ) .* exp( -E0 - E1 * n );