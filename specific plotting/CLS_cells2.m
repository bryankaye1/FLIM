%Plot cells
clear;
%close all;

addpath(fileparts(pwd));
plot_on = 1;

da_cells = [25255,25252,25302,25284,25163,25280]; 
%These segmentations did not look so great: 25395. 25302 was barely
%acceptable.
d_cells = [25265:25270];


representative_samples = [25255,25266,25252,25267];
ind = 0;
for i = da_cells
    ind = ind+1;
    acq_number = [10,5,5,10,10,5];
   % acq_time = [1,1,1,1,1];
    num_pho_thresh=5000;
    matin_wholecell = [28048,28049,28050,28051,28053,28054];

    %pf = 0.0906;
    %Load in data and images
    [input,output,flagoutput] = load_mat_data(i,1,'pause_output_DNE',1,'local',1);
    [spindle] = load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(i)], 'seg_results');
    
    [wholecell_in,wholecell_out,~] = load_mat_data(matin_wholecell(ind),...
        1,'pause_output_DNE',1,'local',1); 
    [whole_cell] = load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(matin_wholecell(ind))], 'seg_results');

  
    %%
    dataname{ind} = output(1,1,1).dataname;
    al = output(1,1,1).w1min/output(1,1,1).w2min;
    %%Get FRET fraction and photons from FRETters
%     for k = 1:length(output)
%     x_group(k) = sum(output(1,1,k).datahis/acq_number(ind));
%     end
%     x_spindle(ind) = sum(x_group);
%    x_cell_simple(ind) = sum(wholecell_in(1,1,1).datahis)/acq_number(ind);
    [x,y,stdpr] = getFF_int(output,al,200000);   
    [x_cell_temp,y_cell,stdpr_cell] = getFF_int(wholecell_out,al,200000);
    
    pixels_cell(ind) = sum(sum(whole_cell.seg_results{4}))/acq_number(ind);
    pixels_spindle(ind) = sum(sum(spindle.seg_results{4}))/acq_number(ind);
    
    
    x = x / ( (pixels_spindle(ind)/8) * acq_number(ind) );
    x_cell(ind) = x_cell_temp / ( pixels_cell(ind) * acq_number(ind) );
    
        figure(ind); clf; subplot(1,2,1); imshow(whole_cell.seg_results{4},...
            'InitialMagnification','fit');
        figure(ind); subplot(1,2,2); imshow(spindle.seg_results{4},...
            'InitialMagnification','fit');
    tub_con = 20;
    F(ind) = wholecell_out.prBest;
    epsilon(ind) = x_cell(ind) / (tub_con*(1-F(ind)+al*F(ind)));   
    [pfm(ind),nmon(ind),fresult] = fit_pf_nm(x,y,stdpr,al);

    if any(i==d_cells)
        pfm(ind) = .0925;
        nmon(ind) = x_int(ind);
    end
    [pol(ind),~,~,polfrac(ind)] = ff2pol(x,y,stdpr,al,pfm(ind),nmon(ind));
    pol_con(ind) = pol(ind)/(epsilon(ind)*length(output));
    %we divide by length of output since pol is a concentration and not a
    %total (you dont sum concentrations to get concentrations)
    
    ci68 = confint(fresult,.682);
    ci95 = confint(fresult,.954);
    
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;
    nmonstd(ind) = (ci68(2,2)-ci68(1,2))/2;
    pfci95(ind) = (ci95(2,1)-ci95(1,1))/2;
    nmonci95(ind) = (ci95(2,2)-ci95(1,2))/2;
    
    
    xmat{ind} = x; %We divide by 8 because we want 
    ymat{ind} = y;
    stdmat{ind} = stdpr;
    
    ti{ind} = strrep(strrep(output(1,1,1).dataname,'.sdt',''),'_',' ');
    thr_pho{ind} = find(xmat{ind}<num_pho_thresh);
    if isfield(output,'blurr')
        ti_im_params{ind} = ['thresh:', num2str(output(1,1,1).im_thr),...
            ' mat:',num2str(i)];
    end
    
    
    pf = pfm(ind);
    polfrac2(ind) = sum( y ./ (pf-y) ) / sum( pf ./ (pf-y ) );
    xmax = max(xmat{ind})*1.3;
    ymax = ymat{ind}(end) + stdmat{ind}(end)+.2;
    if plot_on
        figure(ind); clf; hold on;
        errorbar(xmat{ind},ymat{ind},stdmat{ind},'.');
        errorbar(xmat{ind}(thr_pho{ind}),ymat{ind}(thr_pho{ind}),...
            stdmat{ind}(thr_pho{ind}),'.','Color','red');
        xmodel = 0:xmax*.001:xmax;
        plot(xmodel,pf_nm_pred(al,pfm(ind),nmon(ind),xmodel),'--','Color',...
            [0.6,0.6,0.6]);
        title('FF v int fit','Interpreter','Tex');
        xlabel('intensity');
        ylabel('FRET fraction');
        title(ti{ind});
        axis([0 xmax 0 ymax]);
        figure(2); clf;
        imshow(spindle.seg_results{1}(1:128,1:128),'InitialMagnification','fit');
        title(ti{ind});
    end
end


%[imagedata,seg_results] = load_int_image(round(i));
%image2display = display_segres(imagedata,seg_results); %Delete?
%raw_images{ind} = squeeze(image2display);
%function return pols
%%%%[pol(ind),ave_FRET(ind),ave_FRETerr(ind)] = ff2pol(x_int,y,stdpr,al,pfm(ind));
