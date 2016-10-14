
function [y] = flim_decay_shift(shifti,w2i,backi,sfrc,bins,bneed,irf,Tlaser,wig,data)
intx = 1:sfrc:bins;
irfint = interp1(1:length(irf),irf,intx);
%binskeep = length(data);
binskeep = bins - bneed;

shift = round(shifti/sfrc);
irf2 = circshift(irfint',shift);
if bins<5
ga = 1;
y = ones(length(data),1);
else
ga = interp1(intx, irf2, 1:bins);
s = Tlaser/bins:Tlaser/bins:Tlaser;
f2 = exp(-s/w2i); %signal over one period
f2 = [f2 f2]; %signal over 2 consecutive periods
f2con = conv(f2,ga); %PDF after conv
f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
f2h = f2bar(1:binskeep); %Keep only the appropriate bins
f2h = f2h/sum(f2h);
back = backi/bins;
model = (f2h + back).*wig';
y = model/sum(model)';
end

end