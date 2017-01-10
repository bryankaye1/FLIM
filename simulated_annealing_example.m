%sim_annealing_test
clear;
%Parameter bounds
xmin = 0;
xmax = 100;
Tmin= .001;
Tmax = 1;
T=Tmax;

%cost function
energy = @(x) sin(x)./x;

%starting state
xcurr = rand*xmax +xmin;
xstart= xcurr;
ecurr = energy(xcurr);


ebest = inf;
tsteps = 0;
visited = [];

while T>Tmin
    tsteps = tsteps+1;
    T = T*0.95;
    for i=1:1000
        xpro = xcurr+(rand-0.5)*xmax*sqrt(T/Tmax);
        epro = energy(xpro);
        pacc = exp(-(epro-ecurr)/T);
        if pacc>rand
            xcurr =xpro;
            ecurr = energy(xcurr);
            if ecurr<ebest
                ebest=ecurr;
                xbest = xcurr;
            end
        end
       visited(end+1) = xcurr;
    end
    
    
end
figure(1); clf;
plot(visited,energy(visited)); 
hold on;
plot(xcurr,energy(xcurr),'r.','MarkerSize',10);

