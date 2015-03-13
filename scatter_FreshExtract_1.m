figure(1); clf; figure(2); clf; figure(3); clf; figure(4); clf;

p01 = prBestmat(1,:); pr1 = p01./(1+p01);
p02 = prBestmat(2,:); pr2 = p02./(1+p02); %pr2(18) = mean(pr2(16:17)); pr2(19)=pr2(18);
p03 = prBestmat(3,:); pr3 = p03./(1+p03);

pr1off = pr1(1,:) - mean(pr1(1,end-4:end));
pr2off = pr2(1,:) - mean(pr2(1,end-4:end));
pr3off = pr3(1,:) - mean(pr3(1,end-4:end));

clear x; clear xoff; clear xa; clear xoffa;
for i = 1:3
    
    for j = 1:5
            
            xa(i) = (1.1/100)/2^(i-1); %1.2
            x((i-1)*5+j) = (.8/100)/2^(i-1); %1.2
            xoff((i-1)*5+j) = (.5/100)/2^(i-1); %.75
            xoffa(i) = (.7/100)/2^(i-1); %1.2
    end
end

x = [x,0,0,0,0,0];
xoff = [xoff,0,0,0,0,0];
xa = [xa,0];
xoffa = [xoffa,0];

figure(1); scatter(xoff,pr1off,'b'); hold all; scatter(xoff,pr2off,'g'); scatter(xoff,pr3off,'k'); plot(xoff,xoff,'r');
    title('Fresh Extract Subtracting Offset: 20mg/ml Acc - Buffer Controlled'); 
    xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','+ran later','Slope-1 Line','Location', 'NorthWest');
    
figure(2); scatter(x,pr1(1:end),'b'); hold all;  scatter(x,pr2(1:end),'g'); scatter(x,pr3(1:end),'k'); plot(x,x,'r');
title('Fresh Extract Raw: 20mg/ml Acc - Buffer Controlled'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','+ran later','Slope-1 Line','Location', 'NorthWest');
    
for i =1:4
    pr1a(i) = mean(pr1(5*(i-1)+1:5*(i-1)+5));
    pr2a(i) = mean(pr2(5*(i-1)+1:5*(i-1)+5));
    pr3a(i) = mean(pr3(5*(i-1)+1:5*(i-1)+5));
    
    pr1offa(i) = mean(pr1off(5*(i-1)+1:5*(i-1)+5));
    pr2offa(i) = mean(pr2off(5*(i-1)+1:5*(i-1)+5));
    pr3offa(i) = mean(pr3off(5*(i-1)+1:5*(i-1)+5));
end

figure(3); scatter(xa,pr1a(1:end),'b'); hold all;  scatter(xa,pr2a(1:end),'g'); scatter(xa,pr3a(1:end),'k'); plot(xa,xa,'r');
title('Fresh Extract Raw Averaged: 20mg/ml Acc - Buffer Controlled'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','+ran later','Slope-1 Line','Location', 'NorthWest');
    
figure(4); scatter(xoffa,pr1offa(1:end),'b'); hold all;  scatter(xoffa,pr2offa(1:end),'g'); scatter(xoffa,pr3offa(1:end),'k'); plot(xoffa,xoffa,'r');
title('Fresh Extract Averaged Subtracted Offset: 20mg/ml Acc - Buffer Controlled'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','+ran later','Slope-1 Line','Location', 'NorthWest');
    
    
    