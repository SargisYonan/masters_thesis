%% Map an example of a map for the Stochastic Field Notation Section

figure()

plot(.8, 2.8, 'ro')
hold on

plot(2, 2, 'ro')
hold on

plot(3, 2, 'ro')
hold on

plot(1, 1, 'ro')

xlabel('$x_{1}$','Interpreter','latex')
ylabel('$x_{2}$','Interpreter','latex')
title('Target Field')

axis([.5 3.5 .5 3.5])

set(gca,'xtick',[1:1:1])
set(gca,'ytick',[1:1:1])

grid on
grid minor
