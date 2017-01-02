%Plot cells
clear;
close all;

poster_cells = [25255,25266,25252,25267];
ind = 0;
for i = poster_cells(1)
ind = ind+1;
acq_time = 10;
total_pixels = 128*128;
num_pho_thresh=5000;
pf = 0.0906; 
tic

dataname_ind = {};

%Load in data and images
[~,output,flagoutput] = load_mat_data(i,1,'pause_output_DNE',1);
[imagedata,seg_results] = load_int_image(round(i));

%%
image2display = display_segres(imagedata,seg_results); %Delete?
al = output(1,1,1).w1min/output(1,1,1).w2min;
%%Get FRET fraction and photons from FRETters
[x,y,stdpr] = getFF_int(output,al);
%convert from n_photons to cps (x_int). Fit, fill in missing groups

[x_int,x] = npho2cps(x,total_pixels,acq_time,output);
x_int = x_int/(128*128/8);

[pfm(ind),Nmon(ind)] = fit_pf_nm(x_int,y,stdpr,al);
if isfield(output,'pixel_counts')
    [x_int,x,y,stdpr] = fill_FF_int(x_int,x,y,stdpr,output);
end
xintmat{ind} = x_int;
xmat{ind} = x;
ymat{ind} = y;
stdmat{ind} = stdpr;
raw_images{ind} = squeeze(image2display);
ti{ind} = strrep(strrep(output(1,1,1).dataname,'.sdt',''),'_',' ');
dataname_ind{end+1} = {output(1,1,1).dataname,ind};
thr_pho{ind} = find(xmat{ind}<num_pho_thresh);
if isfield(output,'blurr')
    ti_im_params{ind} = ['thresh:', num2str(output(1,1,1).im_thr),' mat:',num2str(i)];
end
%function return pols
[pol(ind),ave_FRET(ind),ave_FRETerr(ind)] = ff2pol(x_int,y,stdpr,al,pfm(ind));

k = ind;
xmax = 200;
ymax = .1;
figure(1); hold on;
errorbar(xintmat{k},ymat{k},stdmat{k},'.');
errorbar(xintmat{k}(thr_pho{k}),ymat{k}(thr_pho{k}),stdmat{k}(thr_pho{k}),'.','Color','red');
xmodel = 0:xmax*.001:xmax;
plot(xmodel,pf_nm_pred(al,pfm(k),Nmon(k),xmodel),'--','Color',[0.6,0.6,0.6]);
title('FF v int fit','Interpreter','Tex');
xlabel('intensity');
ylabel('FRET fraction');
title(ti{k});
axis([0 xmax 0 ymax]);

figure(2); clf; imshow(seg_results{1}(1:128,1:128)); title(ti{k});

b =1;
end
