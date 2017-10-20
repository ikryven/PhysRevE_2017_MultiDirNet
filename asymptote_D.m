function wa = asymptote_D( u, n )
% wa = asymptote_D( u, n ) computes the asymptote of the size distribution of weakly connected components in
% the directed configuration model.
% u - degree distibution as a two-dimensional array.
% n - one-dimensional array of points where the asymptote is to be calculated.
% Source:  "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
% CC BY 2017. 


excess_distributions;

D = [ 0 -1; 1 0 ];
I = [ 1,   -1   ];
      
s =  -[ m10, m10 ]';

a =      ( ma10 - ma01 ) / 2;
b =  1 - ( ma01 + ma10 ) / 2;


A = ( Sa01 - Sa10 ) / 2;
B = ( Sa01 + Sa10 ) / 2;

cc2 =  a' * D * A * D * a ;
cc1 = (a' * D * B * D * a + a' * D * A * D * b);
cc0 =  a' * D * B * D * b;

if abs( cc2 ) > 1e-10 

    x1 = ( - cc1 + sqrt( cc1^2 - 4 * cc0 * cc2 ) ) / cc2 / 2;
    x2 = ( - cc1 - sqrt( cc1^2 - 4 * cc0 * cc2 ) ) / cc2 / 2;

 if - 1 - eps   < x1 & x1 < 1 + eps
     z = real( x1 );
 elseif -0.5-1  < x2 & x2 < 1 + eps
     z = real( x2 );
 else
     disp('both roots out of the range')
     z=1;
 end;

 if ( -1 < x1 & x1 < 1 ) &  -1 < x2 & x2 < 1
     disp('both roots in the range')
 end

 disp( [x1 x2] );
else
    z = -cc0 / cc1;
    disp( 'only one root' );
    disp( z )
end  

if z < -1
    z = -1;
end

if z > 1
    z = 1;
end;

r_1 = z;

S   =   z * A + B;
aS  = [ S(2,2) -S(1,2); -S(2,1) S(1,1) ];
inS = aS/det(S); 

E_1  = a' * inS * (  a * b' - b * a' ) * inS * b  / ( 2 * a' * inS * a );
E_0  = a' * inS * (  a * b' - b * a' ) * inS * s  / (     a' * inS * a );
E_01 = a' * inS * (  a * s' - s * a' ) * inS * s  / ( 2 * a' * inS * a );

C1 =   ( I * ( Sa10 -  Sa01 ) + ma10' * D * Sa01 - ma01' * D * Sa10 ) * inS;
C2 = - inS * ( Sa01 * D * Sa10 ) * inS;

L0 = 2^(-3/2) *( abs(pi * a'*  aS *a  ) )^(-1/2) ; 
L1 = C1 * s   + (r_1 * a + b )' * ( C2 + C2' ) * s;    
L2 = s' * C2 * s;

wa = L0 * (L1 * n.^(-1.5) + L2 * n.^(-2.5) ) .*  exp( - ( E_1 * n +  E_0 + E_01 ./ n ) );

disp( {'giant component criterion', 2*m10 * m11- m10 * m02  - m10 * m20 + m02 *m20 - m11^2 })



 