figure(2); clf;

plot(TotalInts);
%errorbar(TotalInts(2:end),sqrt(TotalInts(2:end)));

TN1 = TotalInts(2:end)-min(TotalInts(2:end));
TN = TN1/max(TN1);

PN1 = pol(2:end)-min(pol(2:end));
PN = PN1/max(PN1);

figure(3); clf; hold all;
plot(TN(2:end)); plot(PN(2:end));

figure(4); clf; scatter(TN(2:end),PN(2:end));
axis([0,1,0,1]);