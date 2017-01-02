function [y] = pf_nm_pred(al,pf,nmon,x)

fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
y = fitmod(pf,nmon,x);
end