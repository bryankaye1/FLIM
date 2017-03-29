
function [y,varargout] = lifetime_fit(x,bins,bneed,irf,Tlaser,wigs,data,...
    lifetimes,set_w00i)


if lifetimes==2
w2i = x(1);
w02i = x(2);
w1i = x(3);
w01i = x(4);
elseif lifetimes==1
w1i = 0;
w2i = x(1);
w01i = 0;
w02i = x(2); 
else  
fprintf('\n error in number of lifetimes input parameter');
pause;    
end

if ischar(set_w00i)
w00i = 1 - w01i - w02i; 
else
w00i = set_w00i;
end

if w01i + w02i+ w00i > 1
    y = 5e5*(w01i+w02i+w00i-1)+3e3;
else

binskeep = bins - bneed;
binsl = binskeep; %so it loks like cluster script with time removal feature
c = binsl/bins;

s = Tlaser/bins:Tlaser/bins:Tlaser;
f2 = exp(-s/w2i); %signal over one period
f2 = [f2 f2]; %signal over 2 consecutive periods
f2con = conv(f2,irf); %PDF after conv
f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
f2h = f2bar(1:binskeep); %Keep only the appropriate bins
b = sum(f2h)/sum(f2bar);
f2h = f2h/sum(f2h);

f1 = exp(-s/w1i); %signal over one period
f1 = [f1 f1]; %signal over 2 consecutive periods
f1con = conv(f1,irf); %PDF after conv
f1bar = f1con(bins+1:2*bins); %pdf after mod-ing by laser period
f1h = f1bar(1:binskeep); %Keep only the appropriate bins
if sum(f1h)==0
    a = 1;
else    
    a = sum(f1h)/sum(f1bar);
    f1h = f1h/sum(f1h);
end

%back = backi/bins;
%model = (f2h + back).*wig';
%y = model/sum(model)';


rd = 1/(a*w01i+b*w02i+c*w00i);
p2 = data./wigs;

%y = exp(-sum(log((w01i*a*f1h+w02i*b*f2h + w00i/binskeep)*rd).*p2));
y = exp(-sum(log((w01i*f1h+w02i*f2h + w00i/binskeep)*rd).*p2));
%output fitfor plotting etc
varargout{1} = (w01i*f1h+ w02i*f2h + w00i/binskeep)*rd;
varargout{2} = p2;

end
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
