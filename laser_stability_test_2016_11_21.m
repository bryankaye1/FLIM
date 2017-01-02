
clear;
imax = 999;
filename = 'laser_power_test';
pth_sdt = 'C:\Users\Bryan\Documents\MATLAB\data\2016-12-15\';

for i = 1:imax 
    suffix = sprintf('%03s',num2str(i));
    file_name2 = [filename, '_c',suffix,'.sdt'];
    ld = sdt_to_vector(pth_sdt,file_name2);
    dataout(i,:) = ld; 
       
end
%%

total_photons = sum(dataout,2);
acq_time = 50;
dur = 60;
intensity = total_photons / acq_time;

x_plot = (dur:dur:dur*imax)/3600;

% figure; clf; 
% plot(x_plot,intensity,'.','MarkerSize',4);
% xlim([0,x_plot(end)+x_plot(1)/2]);
% xlabel('hours'); ylabel('Counts per second');

int_diff = sum((total_photons(1:end-1,:)-total_photons(2:end,:)).^2,2);
figure(2); clf; 
yyaxis left; plot(x_plot(1:end-1),int_diff);
ylim([0 9e11*.625]);
ylabel('SSE of sequential FLIM decays');
xlabel('hours');
title('Mai-Tai stability at 1000nm');

yyaxis right; plot(x_plot,intensity,'.','MarkerSize',4);
ylabel('total photons');
%ylim([3.8e5 4.6e5]);

for i=2:2:16
figure; clf; hold on;
plot((dataout(60,:)));
plot(dataout(60*i,:)*max(dataout(i,:))/max(dataout(60*i,:)));
xlim([1300,1450]);
ylabel('Counts');
xlabel('Bin Number');
ti = sprintf('IRF hour %1.0f (blue) and %1.0f (orange)',1,i);
title(ti);
end


% plot(intensity,'.','MarkerSize',4);
% xlabel('datapoint'); ylabel('Counts per second');

%120,121