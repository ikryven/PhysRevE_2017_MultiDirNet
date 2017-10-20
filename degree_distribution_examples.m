function u=degree_distribution_examples(id)

if id==2;
    nn=1:20;
    P = 0;
    P( 1, nn ) = exp( -nn *2.266  );
    P( 2, nn ) = exp( -nn * 0.7 );

    P(1,1)=0;
    P = P / sum( P(:) );
    u=P;
    
end



if id==1; 
    [ Y X ] = meshgrid( 0:10, 0:10 );

    c1=[ 0 0 ];
    c2=[ 4 4 ];


    a1=100;
    s1 = 1;
    s2 = 0.4;

    P = a1 * exp( - ( ( X - c1(1) ).^2 + ( Y - c1( 2 ) ).^2 )/s1  )+ exp( - ( ( X - c2(1) ).^2 + ( Y - c2( 2 ) ).^2 )/s2 );


    P(1,1)=0;
    P = P / sum( P(:) );
    u=P;
end;


if id == 3; 
    [ Y X ] = meshgrid( 0:10, 0:10 );

    c1=[ 0 0 ];
    c2=[ 4 4 ];


    a1=220;
    s1 = 0.4;
    s2 = 0.2;

    P = a1 * exp( - ( ( X - c1(1) ).^2 + ( Y - c1( 2 ) ).^2 )/s1  )+ exp( - ( ( X - c2(1) ).^2 + ( Y - c2( 2 ) ).^2 )/s2 );


    P( 1,1 ) = 0;
    P = P / sum( P(:) );
    u = P;
end;
%{

if id==3; 
    [ Y X ] = meshgrid( 0:10, 0:10 );

    c1=[ 1 0 ];
    c2=[ 9 3 ];


    a1=0.9782;

    a2=0.002;

    s1 = 1/5;
    s2 = 1/10;

    P = a1 * exp( - ( ( X - c1(1) ).^2 + ( Y - c1( 2 ) ).^2 )/s1  )+ a2*exp( - ( ( X - c2(1) ).^2 + ( Y - c2( 2 ) ).^2 )/s2 );


    P(1,1)=0;
    P = P / sum( P(:) );

    u=P;

end
%}
