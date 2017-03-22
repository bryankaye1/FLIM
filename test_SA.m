
clear;
for i = 1:1
    a = 5;
localhandle = @(x)SA_testfn(x,a);
options = optimoptions(@simulannealbnd,'FunctionTolerance',0);
[xmin(i),fval(i)] = simulannealbnd(localhandle,90,-100,100,options);%,'FunctionTolerance',{1e-3});
figure(1); clf; hold on;
plot(xmin,fval,'ro');
xplot = linspace(-(max(abs(xmin))+1),(max(abs(xmin(i)))+1),1000);
plot(xplot,localhandle(xplot));
end
mean(fval)
std(fval)