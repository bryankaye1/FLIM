
clear;
imax = 571;
filename = 'irf_40x_air';
pth_sdt = 'C:\Users\Bryan\Documents\MATLAB\data\2016-10-15\';

for i = 1:imax 
    suffix = sprintf('%03s',num2str(i));
    file_name2 = [filename, '_c',suffix,'.sdt'];
    ld = sdt_to_vector(pth_sdt,file_name2);
    dataout(i,:) = ld; 
end
%%
total_photons = sum(dataout,2);
acq_time = 100;
dur = 120;
intensity = total_photons / acq_time;

x_plot = (dur:dur:dur*imax)/60;


figure(1); clf; 
plot(x_plot,intensity,'.','MarkerSize',4);
xlim([0,x_plot(end)+x_plot(1)/2]);
xlabel('minutes'); ylabel('total photon count');
%ylim([0,max(intensity)*1.1]);

figure(2); hold on;
semilogy(dataout(1,:));
semilogy(dataout(end,:));