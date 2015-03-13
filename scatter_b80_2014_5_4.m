
figure(3); clf; figure(4); clf;

p01 = prBestmat(1,1:end); pr1 = p01./(1+p01);
p02 = prBestmat(2,1:end); pr2 = p02./(1+p02);

pr1off = pr1 - mean(pr1(end-3:end));
pr2off = pr2 - mean(pr2(end-3:end));


clear x; clear xoff;
for i = 1:6
    for j = 1:3

            x((i-1)*3+j) = (15.1/100)/2^(i-1);
            xoff((i-1)*3+j) = (14.3/100)/2^(i-1);

    end
end

x = [x,0,0,0];
xoff = [xoff,0,0,0];

figure(3); scatter(xoff,pr1off); hold all; scatter(xoff,pr2off); plot(x,x);
title('BRB80: Subtracting Offset - 20mg/ml Acc (Buffer Controlled)'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
legend('Measurements from Slide 1','Measurements from Slide 2','Slope-1 Line','Location', 'SouthEast');

figure(4); scatter(x,pr1); hold all; scatter(x,pr2); plot(x,x);
title('BRB80: Raw Data - 20mg/ml Acc (Buffer Controlled)'); xlabel('Predicted Acc Fraction'); ylabel('Measured Acc Fraction');
legend('Measurements from Slide 1','Measurements from Slide 2','Slope-1 Line','Location', 'SouthEast');