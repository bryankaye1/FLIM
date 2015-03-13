% 
% di1 = 18;
% di2 = 3;
% 
% for el = 0:di1*di2-1
% fprintf('%f %f \n',ceil((el+1)/di2),mod(el,di2)+1);
% end


a = ones(1,10);
x = 1:1000;
a = x.*x;
a = a/sum(a);
% a2 = fliplr(a);
%  b = cumsum(a2);
%  b2 = fliplr(b)

b2 = fliplr(cumsum(fliplr(a)));

figure(1); clf; hold all;
plot(a,'o');
plot(b2,'o');