% Exapmle2

% In this example, the component size distribution in directed
% configuration model is computed with three different methods:
% fixed point iteration, through convolutions of the degree distribution and
% asymptotically. The asymptote featrures a transient slope -3/2.
% Source:  "Finite connected components in infinite directed and multiplex networks with arbitrary degree distributions " I.Kryven, PhysRevE 2017.
% CC BY 2017.

N = 100;

u = degree_distribution_examples(1);
w1  = fixed_point_directed( u, N*50 ); % for comparision
w  = GL_DirectedNetwork(u,N);
wa = asymptote_D(u,1:N*10);

%%
cla
stairs( w1,  '-r', 'LineWidth',3  );
hold on
stairs( w,  '--k', 'LineWidth',3  );
loglog( wa, '--y', 'LineWidth',3, 'Color', [255, 229, 66 ]/255 );
xlabel('n')
ylabel('w(n)')
ylim([1e-5 1]);
xlim([0 1e3])
set(gca','xscale','log')
set(gca','yscale','log')
legend('fixed-point iteration','GL convolution','asymptote')