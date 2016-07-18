function [new_histogram] = sample_histogram(old_histogram,type,type_param)
%load('C:\Users\Bryan\Documents\MATLAB\FLIM\cdyes\cdata_low_100rep.mat','-mat');
%old_histogram = data(1,1,1).his;
%figure(1); clf;
%plot(old_histogra);%(928:3855));
%set(gca, 'YScale', 'log');

flat_old_histogram = pi*ones(1,sum(old_histogram));
istart = 1;
for i = 1:length(old_histogram)
    flat_old_histogram(istart:istart+old_histogram(i)-1) = i;
    istart = istart + old_histogram(i);
end

if strcmp(type,'fraction')
    n_samples = round(length(flat_old_histogram)*type_param);
elseif strcmp(type,'abs_num')
    n_samples = round(type_param);
else
    warning('did not properly specify sample number, use "fraction" or "abs_num"');
end

edges = 1:length(old_histogram);
indlong = randi(length(flat_old_histogram),1,n_samples); %generate index of photons
samples = flat_old_histogram(indlong); % pull photons
new_histogram = histc(samples, edges); %histogram data
%figure(2);clf; plot(new_histogram);% set(gca, 'YScale', 'log');

end