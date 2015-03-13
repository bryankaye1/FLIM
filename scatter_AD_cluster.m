figure(1); clf; figure(2); clf; figure(3); clf; figure(4); clf;

pr1 = prBest(1,:); %pr1 = p01./(1+p01);
pr2 = prBest(2,:); %pr2 = p02./(1+p02); %pr2(18) = mean(pr2(16:17)); pr2(19)=pr2(18);
%p03 = prBestmat(3,:); pr3 = p03./(1+p03);

pr1off = pr1(1,:) - mean(pr1(1,end-4:end));
pr2off = pr2(1,:) - mean(pr2(1,end-4:end));
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
scatter(xoff,pr1off,'b');
scatter(xoff,pr2off,'g');
title('FE subOff no THR');
xlabel('Rel Acc Con');
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+ Ran','- Ran','Slope-1 Line','Location', 'EastOutside');

figure(2);  hold all; 
scatter(x,pr1(1:end),'b'); 
scatter(x,pr2(1:end),'g');
title('FE no THR'); 
xlabel('Rel Acc Con'); 
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+ Ran','- Ran','Slope-1 Line','Location', 'EastOutside');

for i =1:4
    pr1a(i) = mean(pr1(5*(i-1)+1:5*(i-1)+5));
    pr2a(i) = mean(pr2(5*(i-1)+1:5*(i-1)+5));
    
    pr1offa(i) = mean(pr1off(5*(i-1)+1:5*(i-1)+5));
    pr2offa(i) = mean(pr2off(5*(i-1)+1:5*(i-1)+5));
    
    
    pr1e(i) = std(pr1(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr2e(i) = std(pr2(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    
    pr1offe(i) = std(pr1off(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    pr2offe(i) = std(pr2off(5*(i-1)+1:5*(i-1)+5))/sqrt(5);
    
end

figure(3); hold all;
errorbar(xa,pr1a(1:end),pr1e,':b');
errorbar(xa,pr2a(1:end),pr2e,':r');
title('FE mean no THR');
xlabel('Rel Acc Con');
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+ Ran (SE)' ,'- Ran (SE)','Slope-1 Line','Location', 'EastOutside');

figure(4);hold all;
errorbar(xoffa,pr1offa(1:end),pr1offe,':b');
errorbar(xoffa,pr2offa(1:end),pr2offe,':r');
title('FE mean subOff no THR');
xlabel('Rel Acc Con');
ylabel('Measured Acc Fraction');
axis([0 4 0 max(pr1)]);
legend('+ Ran (SE)','- Ran (SE)','Slope-1 Line','Location', 'EastOutside');


