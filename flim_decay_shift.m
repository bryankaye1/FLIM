
function [y] = flim_decay_shift(x,sfrc,bins,bneed,irf,Tlaser,wigs,data)

shifti = x(1);
w2i = x(2);
w00 = x(3);

intx = 1:sfrc:bins;
irfint = interp1(1:length(irf),irf,intx);
binskeep = bins - bneed;
binsl = binskeep; %so it loks like cluster script with time removal feature
c = binsl/bins;

shift = round(shifti/sfrc);
irf2 = circshift(irfint',shift);

ga = interp1(intx, irf2, 1:bins);
s = Tlaser/bins:Tlaser/bins:Tlaser;
f2 = exp(-s/w2i); %signal over one period
f2 = [f2 f2]; %signal over 2 consecutive periods
f2con = conv(f2,ga); %PDF after conv
f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
f2h = f2bar(1:binskeep); %Keep only the appropriate bins
b = sum(f2h)/sum(f2bar);
f2h = f2h/sum(f2h);

%back = backi/bins;
%model = (f2h + back).*wig';
%y = model/sum(model)';
w02i = 1 - w00;
rd = 1/(b*w02i+c*w00);
p2 = data./wigs';

y = exp(-sum(log((w02i*b*f2h + c*w00/bins)*rd).*p2));

end





% shifti = x(1);
% w2i = x(2);
% backi = x(3);
% 
% intx = 1:sfrc:bins;
% irfint = interp1(1:length(irf),irf,intx);
% %binskeep = length(data);
% binskeep = bins - bneed;
% 
% shift = round(shifti/sfrc);
% irf2 = circshift(irfint',shift);
% if bins<5
% ga = 1;
% y = ones(length(data),1);
% else
% ga = interp1(intx, irf2, 1:bins);
% s = Tlaser/bins:Tlaser/bins:Tlaser;
% f2 = exp(-s/w2i); %signal over one period
% f2 = [f2 f2]; %signal over 2 consecutive periods
% f2con = conv(f2,ga); %PDF after conv
% f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
% f2h = f2bar(1:binskeep); %Keep only the appropriate bins
% f2h = f2h/sum(f2h);
% back = backi/bins;
% model = (f2h + back).*wig';
% y = model/sum(model)';
% end
% 
% %ADD liklihood function here
% exp(sum(log((a*pri*fh + b*w02i*f2h + c*w00/binsl)*rd).*p2));
