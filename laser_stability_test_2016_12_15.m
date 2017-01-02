
clear;
imax = 999;
filename = 'laser_power_test';
pth_sdt = 'C:\Users\Bryan\Documents\MATLAB\data\2016-11-21\';

for i = 1:imax 
    suffix = sprintf('%03s',num2str(i));
    file_name2 = [filename, '_c',suffix,'.sdt'];
    ld = sdt_to_vector(pth_sdt,file_name2);
    irf(i,:) = ld; 
       
end
%%

total_photons = sum(irf,2);
acq_time = 50;
dur = 60;
intensity = total_photons / acq_time;

x_plot = (dur:dur:dur*imax)/3600;

% figure; clf; 
% plot(x_plot,intensity,'.','MarkerSize',4);
% xlim([0,x_plot(end)+x_plot(1)/2]);
% xlabel('hours'); ylabel('Counts per second');

for i = 1:imax-1
    irf_sse(i) = sum( (irf(i,:) - irf(i+1,:)*max(irf(i,:))/max(irf(i+1,:))).^2,2);
end

figure(2); clf; 
yyaxis left; plot(x_plot(1:end-1),irf_sse);
%ylim([0 9e11*.625]);
ylabel('SSE of sequential FLIM decays');
xlabel('hours');
title('Mai-Tai stability at 1000nm');
yyaxis right; plot(x_plot,intensity,'.','MarkerSize',4);
ylabel('total photons');
%ylim([3.8e5 4.6e5]);

figure(1); clf; ind =0;
hour_ref = 1;
for i=hour_ref+1:1:hour_ref+10
subplot(2,5,i-hour_ref); hold on;

plot(1:length(irf(hour_ref*60,:)),1e-5*(irf(hour_ref*60,:)));%,...
  %  '.','MarkerSize',2);
plot(1:length(irf(i*60,:)),1e-5*irf(60*i,:)...
    *max(irf(60*hour_ref,:))/max(irf(60*i,:)));%,'.','MarkerSize',2);
xlim([1275,1400]);
ylabel('Counts (10^{5})');
xlabel('Bin Number');
title('IRF Comparison');
legend(sprintf('IRF hour %1.0f',hour_ref),sprintf('IRF hour %1.0f',i),...
    'Location','south');

end


% plot(intensity,'.','MarkerSize',4);
% xlabel('datapoint'); ylabel('Counts per second');

%120,121