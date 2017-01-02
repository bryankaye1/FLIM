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

sept21_cells_ngr16 = [18634:18649];
sept30 = [20625:20724];
oct6 = 22355:22554;



oct10_donor_nowait = [22556:22605,22656:22705];
oct10_acc_titration = [22606:22655,22706:22755,22756:22805];
oct11 = [22856:22905,22906:22955,22957:23006,23009:23058,23060:23109];
oct11_below_cover = [22806:22855];
oct13 = [23111:23160,23255:23304,23305:23354]-1;

oct13_crashed = [23161:23170,23171:23190,23192:23204,23205:23254]-1;
oct_29 = [23355:23404,23405:23454,23455:23504];

nov1_sp16_5frames = [23505:23515];%,23522:23534];
nov1_sp16_10frames = [23516:23521];

nov1_ngr16_da = [23505:23515];

%fixed = fix_w1(8011:8067,1.8,8238);
ivec = oct11;%oct10_donor_nowait;%[22606:22655,22707:22755,22756: 22805]; %[jul31_donor_ran,jul31_uf1_ran,jul31_uf2_ran];
acq_time = 10;
%%%
show_fit=1;
show_images=1;
tpoint_dur = acq_time+10;
start_time = 120;
total_pixels = 128*128;
num_pho_thresh=5000;
tic
ind = 0;
dataname_ind = {};
for i = ivec
    clear x y stdpr pixel_counts
    ind = ind + 1;
    %Load in data and images
    [~,output,flagoutput] = load_mat_data(i,1);
    imagedata = load_int_image(round(i));
    %%Get FRET fraction and intensity. For now we make x(j)=ni
    ni = output(1,1,1).ni; 
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    for j = 1:length(ni)
        x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
            +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine');
    end
    %convert from n_photons to cps (x_int)
    if isfield(output,'pixel_counts')
        x_int = x.*(total_pixels/acq_time); % x = x.*(128*128/20); % x = ph/pixel  .* pixels/sec
        pixel_counts = output(1,1,1).pixel_counts;
        x = pixel_counts.*x_int.*(acq_time/total_pixels);
    else
        delta_t = acq_time / length(ni);
        x_int = x./delta_t;
    end
    if show_fit && length(ni)>1
        [pfm(ind),Nmon(ind)] = fit_pf_nm(x,y,stdpr,al);
    end
    if isfield(output,'pixel_counts')
        for fill_ind=output(1,1,1).jmax+1:output(1,1,1).num_int_bins
            x_int(fill_ind) = 0;
            x(fill_ind) = 0;
            y(fill_ind) = 0;
            stdpr(fill_ind) = 1e-10;
            pixel_counts(fill_ind)=0;
        end
    end
    xintmat(ind,:) = x_int;
    xmat(ind,:) = x;
    ymat(ind,:) = y;
    stdmat(ind,:) = stdpr;
    raw_images(:,:,ind) = squeeze(imagedata);
    raw_blurred(:,:,ind) = imgaussfilt(imagedata,2);
    ti{ind} = strrep(output(1,1,1).dataname,'.sdt','');
    dataname_ind{end+1} = {output(1,1,1).dataname,ind};
    pf = 0.12;
    
    pol(ind) = sum( x.*y./( pf *( 1 + (al-1).*y) ));
    ave_FRET(ind) = sum(y.*x./(1-y+al*y))/sum(x./(1-y+al*y)); 
    
    yw = (x./(1-y+al*y))/sum(x./(1-y+al*y));
    ywsig = (sqrt(x)./(1-y+al*y))/sum(x./(1-y+al*y));
    FRETerr = y.*yw.*sqrt( (ywsig./yw).^2 + (stdpr./y).^2 );
    ave_FRETerr(ind) = sqrt( sum(FRETerr.^2));
    
    thr_pho{ind} = find(xmat(ind,:)<num_pho_thresh);
end
%Saves images and finds black & white levels
blackval = min(min(min(raw_images(:,:,1:end))));
whiteval = max(max(max(raw_images(:,:,1:end))));
for j = 1:ind
    imseq(:,:,1,j) = mat2gray(raw_images(:,:,j),[blackval whiteval]);
end
%%
%Determine if this is a time series and set appropriate plot parameters
c1index = strfind(strrep(ti{1},'0',''),'_c');
if isempty(c1index)
    x_plot = 1:ind;
    x_limits = [.5,ind+.5];
    xlabel_input = 'sample number';
    subpvec = []
else
    [start_nums,end_nums,~] = find_series(dataname_ind);
    x_plot = [];
    
    %spacer = [60*5,1000,60*8,0]./tpoint_dur;
    spacer = 15;
    ave_FRET1 = (ave_FRET(start_nums(1):end_nums(1)) + ave_FRET(start_nums(2):end_nums(2)))/2;
    ave_FRET2 = (ave_FRET(start_nums(3):end_nums(3)) + ave_FRET(start_nums(4):end_nums(4))...
        +ave_FRET(start_nums(5):end_nums(5)))/3;
    ave_FRET = [ave_FRET1, ave_FRET2];
    
    for i1 = 1:2%length(start_nums)
    x_plot = [x_plot,start_time:tpoint_dur:...
        start_time+tpoint_dur*(end_nums(i1)-start_nums(i1))];
    start_time = start_time+tpoint_dur*(end_nums(i1)-start_nums(i1)+spacer+1);
    end
    x_plot = x_plot./60;
    x_limits = [x_plot(1)-tpoint_dur/2, x_plot(end)+tpoint_dur/2];
    xlabel_input = 'minutes';
end

sumxmat = sum(xmat,2);
F(ind) = struct('cdata',[],'colormap',[]);
fig1 = figure(2);
for j = 1:100
    clf;
    %Place fit in figure (or total photons
    if show_fit && length(ni)>1
        subplot(2,2,1);
        xmax = max(max(xintmat));
        ymax = max(max(ymat));
        hold on;
        errorbar(xintmat(j,:),ymat(j,:),stdmat(j,:),'.');
        errorbar(xintmat(j,thr_pho{j}),ymat(j,thr_pho{j}),stdmat(j,thr_pho{j}),'.','Color','red');
        xmodel = 0:xmax*.001:xmax*1.2;
        plot(xmodel,pf_nm_pred(al,pfm(j),Nmon(j),xmodel),'--','Color',[0.6,0.6,0.6]);
        title('FF v int fit','Interpreter','Tex');
        xlabel('intensity (cps)');
        ylabel('FRET fraction');
        axis([0 xmax*1.2 0 ymax*1.2]);
    else %Or just total number of photons
       % plot(1:j,sumxmat(1:j)','go');        
       % ylabel('Total Photons'); %xlim(x_limits);
    end
    
    %Place intensity image in figure
    pos_im = [.65 .3 .4 .4];
    subplot('Position',pos_im);
    imshow(imseq(:,:,1,j));
    if isempty(c1index)
    title(strrep(ti{j},'_',' '));
    else
    %base_name = strrep(ti{j},'_',' ');
    ctindex = strfind(strrep(ti{j},'0',''),'_c');
    if j<51
        base_name = 'noran';
    else
        base_name = 'with ran';
    end
    %base_name = base_name(1:ctindex);
    tmin = num2str(floor(x_plot(j)));
    tsec = num2str(round(mod(x_plot(j),1)*60),'%02.0f');   
    title([base_name,' ',tmin,':',tsec]);    
    xlabel([num2str(ivec(1)),'-',num2str(ivec(end))])
    end
    
    %Plot total photons, eNpol and mean FRET fraction
    if show_fit && length(ni)>1
    subplot(2,2,3); hold on;
    else
    pos_im = [.1 .1 .5 .8]; subplot('Position',pos_im); hold on;
    end
    yyaxis left;
    plot(x_plot(1:j),sumxmat(1:j)','.','MarkerSize',10);%/(pol(1)/sumxmat(1))
    %plot(x_plot(1:j),pol(1:j),'g.','MarkerSize',10);
    ylim([0,1.1*max(sumxmat)]);
    ylabel('Total Photons');
    
    yyaxis right;
    plot(x_plot(1:j),ave_FRET(1:j),'.','MarkerSize',10);
    ylabel('FRET fraction');
    ylim([0,max(ave_FRET)*1.1]);
%     if max(ave_FRET)>.06
%         ylim([0,max(ave_FRET)*1.1]);
%     end
    xlabel(xlabel_input);
    
    drawnow;
    F(j) = getframe(fig1);
    
end
implay(F);
set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 700 525]);
%movie2avi(F,'C:\Users\Bryan\Documents\MATLAB\movies\oct11_averaged','compression','none');
toc