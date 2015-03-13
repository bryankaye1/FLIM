figure(1); clf; figure(2); clf; figure(3); clf; figure(4); clf;

pr2 = prBest(1,:);
pr4 = prBest(2,:); 
pr1 = prBest(3,:);
pr3 = prBest(4,:); 

pr1off = pr1(1,:) - mean(pr1(1,end-4:end));
pr2off = pr2(1,:) - mean(pr2(1,end-4:end));
pr3off = pr3(1,:) - mean(pr3(1,end-4:end));
pr4off = pr2(1,:) - mean(pr4(1,end-4:end));
%pr3off = pr3(1,:) - mean(pr3(1,end-4:end));

clear x; clear xoff; clear xa; clear xoffa;
for i = 1:3
    
    for j = 1:5
        
        xa(i) = 4/2^(i-1); %1.2
        x((i-1)*5+j) = 4/2^(i-1); %1.2
        xoff((i-1)*5+j) = 4/2^(i-1); %.75
        xoffa(i) = 4/2^(i-1); %1.2
        
    end
end

x = [x,0,0,0,0,0];
xoff = [xoff,0,0,0,0,0];
xa = [xa,0];
xoffa = [xoffa,0];

figure(1); hold all;
scatter(xoff,pr1off);
scatter(xoff,pr2off);
scatter(xoff,pr3off);
scatter(xoff,pr4off);
title('FE subOff');
xlabel('Rel Acc Con');
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+Ran Thr','+Ran/NoThr','-Ran Thr','-Ran/NoThr','Location', 'EastOutside');

figure(2);  hold all; 
scatter(x,pr1(1:end)); 
scatter(x,pr2(1:end));
scatter(x,pr3(1:end)); 
scatter(x,pr4(1:end));
title('FE'); 
xlabel('Rel Acc Con'); 
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+Ran/Thr','+Ran/NoThr','-Ran Thr','-Ran/NoThr','Location', 'EastOutside');

for i =1:4
    pr1a(i) = mean(pr1(5*(i-1)+1:5*(i-1)+5));
    pr2a(i) = mean(pr2(5*(i-1)+1:5*(i-1)+5));
    pr3a(i) = mean(pr3(5*(i-1)+1:5*(i-1)+5));
    pr4a(i) = mean(pr4(5*(i-1)+1:5*(i-1)+5));
    
    pr1offa(i) = mean(pr1off(5*(i-1)+1:5*(i-1)+5));
    pr2offa(i) = mean(pr2off(5*(i-1)+1:5*(i-1)+5));
    pr3offa(i) = mean(pr3off(5*(i-1)+1:5*(i-1)+5));
    pr4offa(i) = mean(pr4off(5*(i-1)+1:5*(i-1)+5));
    
    
    pr1e(i) = std(pr1(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr2e(i) = std(pr2(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr3e(i) = std(pr3(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr4e(i) = std(pr4(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    
    pr1offe(i) = std(pr1off(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr2offe(i) = std(pr2off(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr3offe(i) = std(pr3off(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr4offe(i) = std(pr4off(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    
end

figure(3); hold all;
errorbar(xa,pr1a(1:end),pr1e,':');
errorbar(xa,pr2a(1:end),pr2e,':');
errorbar(xa,pr3a(1:end),pr3e,':');
errorbar(xa,pr4a(1:end),pr4e,':');
title('FE mean');
xlabel('Rel Acc Con');
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+Ran/Thr','+Ran/NoThr','-Ran Thr','-Ran/NoThr','Location', 'EastOutside');

figure(4);hold all;
errorbar(xoffa,pr1offa(1:end),pr1offe,':');
errorbar(xoffa,pr2offa(1:end),pr2offe,':');
errorbar(xoffa,pr3offa(1:end),pr3offe,':');
errorbar(xoffa,pr4offa(1:end),pr4offe,':');
title('FE mean subOff');
xlabel('Rel Acc Con');
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+Ran/Thr','+Ran/NoThr','-Ran Thr','-Ran/NoThr','Location', 'EastOutside');


