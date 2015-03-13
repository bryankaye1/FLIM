

T= 12.58;
tfw = 0.4;
tbac = 0.2;

[~, im] = max(p);

tpb = T/bins;
bfw = round(tfw/tpb);
bbac = round(tbac/tpb);

srem = im - bbac;
erem = im + bfw;

p2 = [p(1:srem),p(erem:end)];

figure(1);clf;plot(p);
hold all; plot(p2);

