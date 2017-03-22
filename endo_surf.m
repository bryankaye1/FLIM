%Plot pf vs unlabeled and acceptor added tubulin to extract.
%3/2017 B.Kaye O.Stiehl


endo = 16;
a = 0:.01:4;
u = 0:.01:45;
m = 2;

[A,U] = meshgrid(a,u);

Z = m * A ./ (A + U + endo);

figure(1); clf;
mesh(A,U,Z);
xlabel('acceptor');
ylabel('unlabeled');

figure(2); clf; hold on;
plot(a,m * (a ./ (a + endo)));
plot(a,m * (a ./ (a +5+ endo)));
plot(a,m * (a ./ (a +10+ endo)));
axis([0 2 0 .25]);

figure(3); clf; hold on;
plot(u,m * (2 ./ (2 +u+ endo)));
plot(u,m * (1 ./ (1 +u+ endo)));
xticks(0:7.5:45);