
x = 10:10:1000;
p01 = prBestmat(1,:); pr1 = p01./(1+p01);
p02 = prBestmat(2,:); pr2 = p02./(1+p02);


figure(3); clf; hold all;
plot(x,p01,'b');  
plot(x,p02,'g'); 
title('Time Series');
xlabel('seconds'); ylabel('FRET fraction');
legend('TS1 1.5e5 ph/s','TS2 4.5e5 ph/s');


