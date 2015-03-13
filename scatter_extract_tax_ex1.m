figure(1); clf; figure(2); clf;

p01 = prBestmat(1,:); pr1 = p01./(1+p01);
p02 = prBestmat(2,:); pr2 = p02./(1+p02);

pr1off = pr1(1,:) - mean(pr1(1,end-2:end));
pr2off = pr2(1,:) - mean(pr2(1,end-2:end));


clear x; clear xoff;
for i = 1:4
    for j = 1:3

            x((i-1)*3+j) = (5.5/100)/2^(i-1);
            xoff((i-1)*3+j) = (5/100)/2^(i-1);

    end
end

x = [x,0,0,0];
xoff = [xoff,0,0,0];
figure(1); scatter(xoff,pr1off); hold all; scatter(xoff,pr2off); plot(xoff,xoff);
title('Extract-Taxol: Subtracting Offset - Acc 20mg/ml (buffer controlled)');
xlabel('Predicted Acc Fraction'); ylabel('Measured Acceptor Fraction');
legend('Measurement on Slide 1','Measurement on Slide 2','Slope-1 Line','Location', 'SouthEast');

figure(2); scatter(x,pr1(1:end)); hold all; scatter(x,pr2(1:end)); plot(x,x);
title('Extract-Taxol: Raw Data - Acc 20mg/ml (buffer controlled)');
xlabel('Predicted Acc Fraction'); ylabel('Measured Acceptor Fraction');
legend('Measurement on Slide 1','Measurement on Slide 2','Slope-1 Line','Location', 'SouthEast');