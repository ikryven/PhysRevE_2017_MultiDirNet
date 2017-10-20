function [mu S] = make_S_mu( u );

[ Y X ] = meshgrid( 0:(size(u,2)-1),0:(size(u,1)-1) );


mu = [sum( sum( X .* u )); sum( sum( Y .* u ))];

s_11 = sum( sum( u .* ( X - mu( 1 ) ).^2 ));
s_22 = sum( sum( u .* ( Y - mu( 2 ) ).^2 ));

s_12 = sum( sum( u .* ( X - mu( 1 ) ) .* ( Y - mu( 2 ) ) ) );
s_21 = sum( sum( u .* ( Y - mu( 2 ) ) .* ( X - mu( 1 ) ) ) );

S = [s_11 s_12; s_21 s_22];