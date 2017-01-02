%Calculate Polymer:

acc = 250;
acclab = acc*.66;
donor = 250;
ex = 15;

acv = [3,2.5,2,1.5,1,0.5];

for i = 1:6
    m(i) = acclab*acv(i) / (acc*acv(i)+1*donor+40*ex);
end
pf = 2.*m'