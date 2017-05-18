%Plot cells
clear;
%close all;

addpath(fileparts(pwd));
plot_on = 1;

da_cells = [25255,25252,25302,25284,25395,25163,25280];
d_cells = [25265:25270];


representative_samples = [25255,25266,25252,25267];
ind = 0;
for i = da_cells
    ind = ind+1;
    acq_time = [10,5,5,10,5,10,5];
    total_pixels = 128*128;
    num_pho_thresh=5000;
    matin_wholecell = [28034:28040];

    %pf = 0.0906;
    %Load in data and images
    [input,output,flagoutput] = load_mat_data(i,1,'pause_output_DNE',1,'local',1);
    [seg_results_cell] = load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(matin_wholecell(ind))], 'seg_results');
    [whole_cell] = load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(matin_wholecell(ind))], 'input');
    [seg_results_spindle] = load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(i)], 'seg_results');
  
    %%
    dataname{ind} = output(1,1,1).dataname;
    al = output(1,1,1).w1min/output(1,1,1).w2min;
    %%Get FRET fraction and photons from FRETters
    [x,y,stdpr] = getFF_int(output,al,200000);
    %convert from n_photons to cps (x_int). Fit, fill in missing groups
    
    [x_int,x] = npho2cps(x,total_pixels,acq_time(ind),output);
    x_int = x_int/(128*128/8);
    
    cell_pho = whole_cell.input(1,1,1).ni;
    cell_int = npho2cps(cell_pho,total_pixels,acq_time(ind),output);
    cell_int = cell_int/(128*128); %divide by 8?
    
    pixels_cell = sum(sum(seg_results_cell.seg_results{4}));
    pixels_spindle = sum(sum(seg_results_spindle.seg_results{4}));
    
    tub_con = 20;
    v_cell = pixels_cell/pixels_spindle;
    epsilon = cell_int / (tub_con * v_cell);
    
    [pfm(ind),nmon(ind),fresult] = fit_pf_nm(x_int,y,stdpr,al);
    ci68 = confint(fresult,.682);
    ci95 = confint(fresult,.954);
    
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;
    nmonstd(ind) = (ci68(2,2)-ci68(1,2))/2;
    pfci95(ind) = (ci95(2,1)-ci95(1,1))/2;
    nmonci95(ind) = (ci95(2,2)-ci95(1,2))/2;
    
    
    xintmat{ind} = x_int;
    xmat{ind} = x;
    ymat{ind} = y;
    stdmat{ind} = stdpr;
    
    ti{ind} = strrep(strrep(output(1,1,1).dataname,'.sdt',''),'_',' ');
    thr_pho{ind} = find(xmat{ind}<num_pho_thresh);
    if isfield(output,'blurr')
        ti_im_params{ind} = ['thresh:', num2str(output(1,1,1).im_thr),...
            ' mat:',num2str(i)];
    end
    
    if any(i==d_cells)
        pfm(ind) = .0925;
        nmon(ind) = x_int(ind);
    end
    [pol(ind),~,~,polfrac(ind)] = ff2pol(x_int,y,stdpr,al,pfm(ind),nmon(ind));
    pol_ep(ind) = pol(ind)/epsilon;
    
    pf = pfm(ind);
    polfrac2(ind) = sum( y ./ (pf-y) ) / sum( pf ./ (pf-y ) );
    xmax = max(xintmat{ind})*1.3;
    ymax = ymat{ind}(end) + stdmat{ind}(end)+.2;
    if plot_on
        figure(ind); clf; hold on;
        errorbar(xintmat{ind},ymat{ind},stdmat{ind},'.');
        errorbar(xintmat{ind}(thr_pho{ind}),ymat{ind}(thr_pho{ind}),...
            stdmat{ind}(thr_pho{ind}),'.','Color','red');
        xmodel = 0:xmax*.001:xmax;
        plot(xmodel,pf_nm_pred(al,pfm(ind),nmon(ind),xmodel),'--','Color',...
            [0.6,0.6,0.6]);
        title('FF v int fit','Interpreter','Tex');
        xlabel('intensity');
        ylabel('FRET fraction');
        title(ti{ind});
        
        
    end
    %axis([0 xmax 0 ymax]);
    
    %figure(2); clf; imshow(seg_results{1}(1:128,1:128)); title(ti{k});
end


%[imagedata,seg_results] = load_int_image(round(i));
%image2display = display_segres(imagedata,seg_results); %Delete?
%raw_images{ind} = squeeze(image2display);
%function return pols
%%%%[pol(ind),ave_FRET(ind),ave_FRETerr(ind)] = ff2pol(x_int,y,stdpr,al,pfm(ind));
