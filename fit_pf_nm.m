function [pfm,nmon] = fit_pf_nm(x,y,stdpr,al)

if min(x)==0
istart = find(x==0,1,'last')+1;
x = x(istart:end);
y = y(istart:end);
stdpr = stdpr(istart:end);
end
fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
fresult = fit((x)',y',fitmod,'StartPoint',[max(y),min(x)],'Weights',...
    (1./(stdpr)).^2','Lower',[0,0],'Upper',[1,max(x)]);

pfm = fresult.a;
nmon = fresult.b;

end