figure(1); clf; figure(2); clf; figure(3); clf; figure(4); clf;

imax =3;
jmax =5;

p01 = prBestmat(1,:,2); pr1 = p01./(1+p01); %pr1(2) = pr1(1);
p02 = prBestmat(2,:,2); pr2 = p02./(1+p02);
% p03 = prBestmat(1,:,2); pr3 = p03./(1+p03);
% p04 = prBestmat(2,:,2); pr4 = p04./(1+p04);

%p03 = prBestmat(3,:); pr3 = p03./(1+p03);

pr1off = pr1(1,:) - mean(pr1(1,end-jmax+1:end));
pr2off = pr2(1,:) - mean(pr2(1,end-jmax+1:end));
%pr3off = pr3(1,:) - mean(pr3(1,end-4:end));

clear x; clear xoff; clear xa; clear xoffa;


for i = 1:3
    
    for j = 1:jmax
            
            xa(i) = (1.6/100)/2^(i-1); %1.2
            x((i-1)*jmax+j) = (1.6/100)/2^(i-1); %1.2
            xoff((i-1)*jmax+j) = (1.1/100)/2^(i-1); %.75
            xoffa(i) = (1.2/100)/2^(i-1); %1.2
    end
end

x = [x,zeros(1,jmax)];
xoff = [xoff,zeros(1,jmax)];
xa = [xa,0];
xoffa = [xoffa,0];

figure(1); scatter(xoff,pr1off,'b'); hold all; scatter(xoff,pr2off,'g'); plot(xoff,xoff,'r');
    title('Fresh Extract Subtracting Offset: 20mg/ml Acc - Buffer Controlled'); 
    xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','Slope-1 Line','Location', 'NorthWest');
    
figure(2); scatter(x,pr1(1:end),'b'); hold all;  scatter(x,pr2(1:end),'g');  plot(x,x,'r');
title('Fresh Extract Raw: 20mg/ml Acc - Buffer Controlled'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','Slope-1 Line','Location', 'NorthWest');
    
for i =1:imax+1
    pr1a(i) = mean(pr1(jmax*(i-1)+1:jmax*(i-1)+jmax));
    pr2a(i) = mean(pr2(jmax*(i-1)+1:jmax*(i-1)+jmax));
   
    
    pr1offa(i) = mean(pr1off(jmax*(i-1)+1:jmax*(i-1)+jmax));
    pr2offa(i) = mean(pr2off(jmax*(i-1)+1:jmax*(i-1)+jmax));

end

figure(3); scatter(xa,pr1a(1:end),'b'); hold all;  scatter(xa,pr2a(1:end),'g');  plot(xa,xa,'r');
title('Fresh Extract Raw Averaged: 20mg/ml Acc - Buffer Controlled'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','Slope-1 Line','Location', 'NorthWest');
    
figure(4); scatter(xoffa,pr1offa(1:end),'b'); hold all;  scatter(xoffa,pr2offa(1:end),'g'); plot(xoffa,xoffa,'r');
title('Fresh Extract Averaged Thresholded Subtracted Offset: 20mg/ml Acc - Buffer Controlled'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
    legend('+ Ran','- Ran','Slope-1 Line','Location', 'NorthWest');
    
    
    