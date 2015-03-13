
x = 70:10:560;
p01 = prBestmat(1,:,1); pr1 = p01./(1+p01);
p02 = prBestmat(1,:,2); pr2 = p02./(1+p02);


figure(3); clf; hold all;
plot(x,p01,'b');  
%plot(x,p02,'g'); 
title('Time Series of Ran Aster Nucleation');
xlabel('seconds'); ylabel('FRET Fraction');
% legend('TS2 ~4e5: ph/s','TS1 ~1e4: ph/s','Location','Southwest');
%set(gca,'XTick',[70,200:200:800,1070]);

smp01 = smooth(p01,10);
smp02 = smooth(p02,10);
figure(4); clf; hold all;
plot(x,smp01,'b');  
%plot(x,smp02,'g'); 
title('Time Series Smoothed over 10 timepoints');
xlabel('seconds'); ylabel('FRET Fraction');
% legend('TS2 ~4e5: ph/s','TS1 ~1e4: ph/s','Location','Southwest');
%set(gca,'XTick',[70,200:200:800,1070]);

% fprintf('\nTS1 Z-score : %f (STD/MEAN)\n', std(p01)/mean(p01));
% fprintf('\nTS2 Z-score : %f (STD/MEAN)\n', std(p02)/mean(p02));
% fprintf('\nZ-TS1 / Z-TS2 : %f , sqrt(3) is %f\n', std(p02)/mean(p02),1/sqrt(3));