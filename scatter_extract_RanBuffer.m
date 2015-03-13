figure(1); clf; figure(2); clf;

p01 = prBestmat(1,:); pr1 = p01./(1+p01);
p02 = prBestmat(2,:); pr2 = p02./(1+p02);

pr1off = pr1(1,:) - mean(pr1(1,end-4:end));
pr2off = pr2(1,:) - mean(pr2(1,end-4:end));


clear x; clear xoff;
for i = 1:3
    for j = 1:5

            x((i-1)*5+j) = (1.2/100)/2^(i-1);
            xoff((i-1)*5+j) = (.8/100)/2^(i-1);

    end
end

x = [x,0,0,0,0,0];
xoff = [xoff,0,0,0,0,0];
figure(1); scatter(xoff,pr1off,'b'); hold all;  plot(xoff,xoff); scatter(xoff,pr2off,'g');
title('Extract-Ran Buffer Subtracting Offset: 20mg/ml (Buffer Controlled)'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
legend('irf-before','Slope-1 Line','irf-after','Location', 'SouthEast');

figure(2); scatter(x,pr1(1:end),'b'); hold all;  plot(x,x); scatter(x,pr2(1:end),'g');
title('Extract-Ran Raw: 20mg/ml (Buffer Controlled)'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
legend('irf-before','Slope-1 Line','irf-after','Location', 'SouthEast');