clear;
m3 = load('m3data');
m2 = load('m2data');
m1 = load('m1data');
m0 = load('m0data');

m0norm3 = sum(m3.p(2700:end))/sum(m0.p(2700:end));
m0n3 = m0.p*m0norm3;

m0norm2 = sum(m2.p(2700:end))/sum(m0.p(2700:end));
m0n2 = m0.p*m0norm2;

m0norm1 = sum(m1.p(2700:end))/sum(m0.p(2700:end));
m0n1 = m0.p*m0norm1;

figure(1);clf;hold on;
plot(log(m0n1));
plot(log(m1.p));

figure(2);clf;hold on;
plot(log(m0n2));
plot(log(m2.p));

figure(3);clf;hold on;
plot(log(m0n3));
plot(log(m3.p));

% figure(3);
% plot(log(m0n));
% hold all;
% plot(log(m3.p));