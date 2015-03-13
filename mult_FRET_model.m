
for i = 0:2
s = 3;
a = .2/2^i;

s2 = s-1;
na = 1-a;

%pg2 = 1-(1-a)^s-s*a*(1-a)^(s-1)
%pf2 = 1-(na+s*a)*(na)^(s-1)
pg2 = 1 - (a*s2+1)*(1-a)^s2;


prob(1+i) = a*s*(1-a)^s2-pg2;

end
prob