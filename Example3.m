% Exapmle2

% In this example, the component size distribution in directed
% configuration model is computed with three different methods:
% fixed point iteration, through convolutions of the degree distribution, and
% asymptotically. The asymptote featrures a transient slope -1/2.
% Source:  "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
% CC BY 2017.

N = 100;

u = degree_distribution_examples( 2 );
w1 = component_num_degenerate( u, N*100 );
w =  components_lagrange_DD(u,N*10);
wa = asymptote_DD( u, 1:N * 10 );

%%
cla
stairs( w1,  '-r', 'LineWidth',3  );
stairs( w,  '--k', 'LineWidth',3  );
hold on
loglog( wa, '--y', 'LineWidth',3, 'Color', [255, 229, 66 ]/255 );

xlabel( 'n'    );
ylabel( 'w(n)' );

ylim([1e-4 2e-2]);
xlim([0 1e3])
set(gca','xscale','log')
set(gca','yscale','log')
legend('fixed-point iteration','GL convolution','asymptote')