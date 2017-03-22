


shifti = rand*shiftmax+shiftmin;
w2i = rand*w2max+w2main;

backi = backmin:backstep:backmax


shift = round(shifti/sfrac);
irf2 = circshift(irfint',shift);
ga = interp1(intx, irf2, 1:bins);


s = Tlaser/bins:Tlaser/bins:Tlaser;
f2 = exp(-s/w2i); %signal over one period
f2 = [f2 f2]; %signal over 2 consecutive periods
f2con = conv(f2,ga); %PDF after conv
f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
f2h = f2bar(1:binskeep); %Keep only the appropriate bins
f2h = f2h/sum(f2h);


back = backi/bins;
model = (f2h + back).*wigs;

model = model/sum(model);
data = data1/sum(data1);
ga =ga/sum(ga);

res = (data-model)./sqrt(model); %changed from sqrt(data) to sqrt(model)
sumresi = sum(res.^2);

if sumresi < sumres
    sumres = sumresi;
    shiftb = shifti;
    gab = ga;
    w2b = w2i;
    backb = backi;
    wigshiftb = wigshifti;
    wigsb = wigs;
    
    modelb = model;
    datab = data;
    resb = res;
    