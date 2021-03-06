%%change nph_counts to pixel_counts after I rerun the int_binned analysis
%In the fit function, I added in a regularization to the weighting to
%(heuristically) account for xvalue uncertainty

%close all;
clear;
%june 23 w1sweep @ 2e5 sp=16:  [8931:8987,8988:9044,9045:9101,9102:9158];
%june 23 w1sweep @ 2e5 ngr=16:  [8011:8067,8068:8124,8125:8181,8182:8238]

aug30_uf1_noran = [16384:16443,16323:16382];%,16383]; %noran
aug30_uf2_ran = [16504:16564,16444:16503]; %ran

jul31_donor_ran = 15109:15138;
jul31_uf2_ran = 12549:12608;
jul31_uf1_ran = 12507:12548;

july31_UF_endpoints = 15057:15066;
jul31_dr_endpoints = 15095:15102;

aug19_UF4_ran= 15549:15608;
aug19_UF4_noran= 15489:15548;
aug19_UF3_noran= 15369:15428;
aug19_UF2_ran= 15224:15248;
aug19_UF1_noran= 15199:15223;

sept6_2e5 = [17120:17269,17494:17543];
sept9 = 17544:17692;
sept13 = [18092:18141,18042:18091,18142:18191,17992:18041,17841:17891,17942:17991,17892:17941];
sept14_old_irf = [18192:18241,18242:18291,18292:18341];
sept14 = [18342:18391,18492:18541,18392:18441,18542:18591,18442:18491];
sept30 = [20625:20724];
oct6 = 22355:22554;

sept21_cells_ngr16 = [18634:18649];

oct10_donor_nowait = [22556:22605,22656:22705];
oct10_acc_titration = [22606:22655,22706:22755,22756:22805];
oct11 = [22856:22905,22906:22955,22957:23006,23009:23058,23060:23109];
oct11_below_cover = [22806:22855];
oct13 = [23111:23160,23255:23304,23305:23354]-1;

oct13_crashed = [23161:23170,23171:23190,23192:23204,23205:23254]-1;
oct_29 = [23355:23404,23405:23454,23455:23504];

nov1_sp16_5frames = [23505:23515];%,23522:23534];
nov1_sp16_10frames = [23516:23521];

nov1_ngr16_da = [23535:23545];
nov1_bestcells = [23614:23627]; %segmented cells nsp 8, blurr =1, thresh = 0.25
cherry_cells = [25225,25255,25315,25345,25252,25312,25302,25332];
reg_cells = [25284,25314,25395,25163,25280];
poster_cells = [25255,25252,25266,25267];

fig5_control = [25078:25090,25783:25795];
fig3_control = [25079,25784,25080,25785,25078,25783];

%fixed = fix_w1(8011:8067,1.8,8238);
ivec = fig3_control;%[22606:22655,22707:22755,22756: 22805]; %[jul31_donor_ran,jul31_uf1_ran,jul31_uf2_ran];
acq_time = 10;
%%%
show_fit=1;
show_images=1;
tpoint_dur = acq_time+10;
start_time = 120;
total_pixels = 128*128;
num_pho_thresh=5000;
pf = 0.12;
tic
ind = 0;
dataname_ind = {};

for i = ivec
    clear x y stdpr pixel_counts
    ind = ind + 1;
    %Load in data and images
    [~,output,flagoutput] = load_mat_data(i,1,'pause_output_DNE',1);
    [imagedata,seg_results] = load_int_image(round(i));
    image2display = display_segres(imagedata,seg_results); %Delete?
    al = output(1,1,1).w1min/output(1,1,1).w2min;
    %%Get FRET fraction and photons from FRETters
    [x,y,stdpr] = getFF_int(output,al);
    %convert from n_photons to cps (x_int). Fit, fill in missing groups

    [x_int,x] = npho2cps(x,total_pixels,acq_time,output);
    if show_fit && length(x)>1
        [pfm(ind),Nmon(ind)] = fit_pf_nm(x_int,y,stdpr,al);
    end 
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
    [pol(ind),ave_FRET(ind),ave_FRETerr(ind)] = ff2pol(x,y,stdpr,al,pf);
    
end


%%
%Determine if this is a time series and set appropriate plot parameters
[start_nums,end_nums,~] = find_series(dataname_ind);
if isempty(start_nums)
    x_plot = 1:ind;
    x_limits = [.5,ind+.5];
    xlabel_input = 'sample number';
else
    x_plot = [];
    %spacer = [60*5,1000,60*8,0]./tpoint_dur;
    spacer = 15;
    for i1 = 1:length(start_nums)
    x_plot = [x_plot,start_time:tpoint_dur:...
        start_time+tpoint_dur*(end_nums(i1)-start_nums(i1))];
    start_time = start_time+tpoint_dur*(end_nums(i1)-start_nums(i1)+spacer+1);
    end
    x_plot = x_plot./60;
    x_limits = [x_plot(1)-tpoint_dur/2, x_plot(end)+tpoint_dur/2];
    xlabel_input = 'minutes';
end

%sumxmat = sum(xmat,2);
F(ind) = struct('cdata',[],'colormap',[]);
fig1 = figure(2);
for k = 1:ind
   % if Nmon(k)>2e4
    clf;
    %Place fit in figure (or total photons
    if show_fit && length(x)>1
        subplot(1,2,1);
        xmax = max(xintmat{k});
        hold on;
        errorbar(xintmat{k},ymat{k},stdmat{k},'.');
        errorbar(xintmat{k}(thr_pho{k}),ymat{k}(thr_pho{k}),stdmat{k}(thr_pho{k}),'.','Color','red');
        xmodel = 0:xmax*.001:xmax*1.2;
        plot(xmodel,pf_nm_pred(al,pfm(k),Nmon(k),xmodel),'--','Color',[0.6,0.6,0.6]);
        title('FF v int fit','Interpreter','Tex');
        xlabel('intensity (cps)');
        ylabel('FRET fraction');
        title(ti{k});
        %axis([0 xmax*1.2 0 ymax*1.2]);
    else %Or just total number of photons
       % plot(1:j,sumxmat(1:j)','go');        
       % ylabel('Total Photons'); %xlim(x_limits);
    end
    
    %Place intensity image in figure
    pos_im = [.6 .3 .4 .4];
    %subplot('Position',pos_im);
    subplot(1,2,2);
    imshow(raw_images{k});
    if isempty(start_nums)
    title(strrep(ti{k},'_',' '));
    else
    base_name = strrep(ti{k},'_',' ');
    ctindex = strfind(strrep(ti{k},'0',''),'_c');
    base_name = base_name(1:ctindex);
    tmin = num2str(floor(x_plot(k)));
    tsec = num2str(round(mod(x_plot(k),1)*60),'%02.0f');   
    title([base_name,' ',tmin,':',tsec]);    
    xlabel([num2str(ivec(1)),'-',num2str(ivec(end))])
    end
    
    %Plot total photons, eNpol and mean FRET fraction
%     if show_fit && length(x)>1
%     subplot(2,2,3); hold on;
%     else
%     pos_im = [.1 .1 .5 .8]; subplot('Position',pos_im); hold on;
%     end
%     yyaxis left;
%     plot(x_plot(1:j),sumxmat(1:j)','.','MarkerSize',10);%/(pol(1)/sumxmat(1))
%     %plot(x_plot(1:j),pol(1:j),'g.','MarkerSize',10);
%     ylim([0,1.1*max(sumxmat)]);
%     ylabel('Total Photons');
%     
%     yyaxis right;
%     plot(x_plot(1:j),ave_FRET(1:j),'.','MarkerSize',10);
%     ylabel('FRET fraction');
%     ylim([0,max(ave_FRET)*1.1]);
% %     if max(ave_FRET)>.06
% %         ylim([0,max(ave_FRET)*1.1]);
% %     end
%     xlabel(xlabel_input);
    
    drawnow;
    F(k) = getframe(fig1);
    %end
end
implay(F);
set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 1000 800]);


%movie2avi(F,'/Users/bryankaye/Desktop/thaw','compression','none');
toc