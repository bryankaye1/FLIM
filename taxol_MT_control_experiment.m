e1 = [0.009093333, 0.016242222, 0.023064444,  0.034262222, 0.054611111, 0.091034444];
e2 = [0.007264444, 0.019617778, 0.02278, 0.035328889, 0.055287778, 0.093141667];
ea = e1 + e2;

pred = [0, 1.25, 2.5, 5, 10, 20];
pred = pred*0.63/100;

% figure(13); clf;
% plot(log(pred), log(ea));


p0norm = p0/max(p0);
p1norm = p1/max(p1);
p2norm = p2/max(p2);
p3norm = p3/max(p3);
p4norm = p4/max(p4);
p5norm = p5/max(p5);


p0n = (p0-min(p0));
p5n = (p5-min(p5));



tkeepi = tmaxi-tmini;
Tfrac = tkeepi/4096;
laser_period = 12.58;
recwin = laser_period*Tfrac;

tpb = recwin/length(p0);
time = tpb:tpb:recwin;



figure(1); clf; plot(time,log(p0norm)); hold on;
% figure(2); clf; plot(time,log(p1norm)); pause(0.3);
% figure(3); clf; plot(time,log(p2norm)); pause(0.3);
% figure(4); clf; plot(time,log(p3norm)); pause(0.3);
% figure(5); clf; plot(time,log(p4norm)); pause(0.3);
plot(time,log(p5norm)); pause(0.3);




% 
% p0m = (p0-p0(end));
% p5m = (p5-p5(end));


% 
% figure(2); 
% clf; plot(time,log(p0n));  hold on;
% plot(time,log(p5n), 'r'); 
% 
% 
% figure(3);
% clf; plot(time(1:tkeepi-713),log(p0m(715:end)));  hold on;
% plot(time(1:tkeepi-713),log(p5m(715:end)), 'r'); 



