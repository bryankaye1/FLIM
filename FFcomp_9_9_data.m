%%change nph_counts to pixel_counts after I rerun the int_binned analysis
%In the fit function, I added in a regularization to the weighting to
%(heuristically) account for xvalue uncertainty

%close all;
clear;
aug30_uf1_noran = [16384:16443,16323:16382];%,16383]; %noran
aug30_uf2_ran = [16504:16564,16444:16503]; %ran

jul31_donor_ran = 15109:15138;
jul31_uf2_ran = 12549:12608;
jul31_uf1_ran = 12507:12548;

jul31_dr_endpoints = 15095:15102;

aug19_UF4_ran= 15549:15608;
aug19_UF4_noran= 15489:15548;
aug19_UF3_noran= 15369:15428;
aug19_UF2_ran= 15224:15248;
aug19_UF1_noran= 15199:15223;

sept6_2e5 = 17120:17269;
sept9 = 17544:17692;
sept9_noranshift = 17693:17792;

ivec = sept9_noranshift; %[jul31_donor_ran,jul31_uf1_ran,jul31_uf2_ran];
acq_time = 20;
%%%
show_fit=1;
show_images=1;
save_time = 30;
start_time = 45;
total_pixels = 128*128;
num_pho_thresh=5000;
tic
ind = 0;
for i = ivec
    clear x y stdpr pixel_counts
    ind = ind + 1;
    %Load in data and images
    [~,output,~] = load_mat_data(i,1);
    imagedata = load_int_image(round(i));
    %%Determine if Super Pixels (pixc) are present
    %     try pixel_counts = output(1,1,1).pixel_counts;
    %         pixc = 1;
    %     catch
    %         pixc = 0;
    %     end
    %%Get FRET fraction and intensity. For now we make x(j)=ni
    ni = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    for j = 1:length(ni);
        x(j) = ni(j);%*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
        %+sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'combine_background');
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
    subpvec = [];
else
    numsam1 = length([17544:17593]);
    numsam2 = length([17594:17643]);
    numsam3 = length([17644:17692]);
    x_plot1 = (start_time:save_time:start_time+save_time*(numsam1-1))./60;
    x_plot2 = (start_time:save_time:start_time+save_time*(numsam2-1))./60;
    x_plot3 = (start_time:save_time:start_time+save_time*(numsam3-1))./60;
    x_plot = [x_plot1,x_plot2,x_plot3];
    x_limits = [x_plot(1)-save_time/2, x_plot(end)+save_time/2];
    xlabel_input = 'minutes';
end

sumxmat = sum(xmat,2);
F(ind) = struct('cdata',[],'colormap',[]);
fig1 = figure(2); clf;
for j = 1:ind
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
        base_name = strrep(ti{j},'_',' ');
        ctindex = strfind(strrep(ti{j},'0',''),'_c');
        base_name = base_name(1:ctindex);
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
    if j<=numsam1
        plot(x_plot(1:j),sumxmat(1:j),'b.','MarkerSize',10);
    elseif j<=(numsam2+numsam1)
        plot(x_plot(1:numsam1),sumxmat(1:numsam1),'b.','MarkerSize',10);
        plot(x_plot(numsam1+1:j),sumxmat(numsam1+1:j),'g.','MarkerSize',10);
    else
        plot(x_plot(1:numsam1),sumxmat(1:numsam1),'b.','MarkerSize',10);
        plot(x_plot(numsam1+1:j),sumxmat(numsam1+1:j),'g.','MarkerSize',10);
        plot(x_plot(numsam2+numsam1+1:j),sumxmat(numsam1+numsam2+1:j),'r.','MarkerSize',10);
    end
    ylim([0,1.05*max(sumxmat)]);
    ylabel('Total Photons');
    
    yyaxis right;
    if j<=numsam1
        plot(x_plot(1:j),ave_FRET(1:j),'b.','MarkerSize',10);
    elseif j<=(numsam2+numsam1)
        plot(x_plot(1:numsam1),ave_FRET(1:numsam1),'b.','MarkerSize',10);
        plot(x_plot(numsam1+1:j),ave_FRET(numsam1+1:j),'g.','MarkerSize',10);
    else
        plot(x_plot(1:numsam1),ave_FRET(1:numsam1),'b.','MarkerSize',10);
        plot(x_plot(numsam1+1:j),ave_FRET(numsam1+1:j),'g.','MarkerSize',10);
        plot(x_plot(numsam2+numsam1+1:j),ave_FRET(numsam2+numsam1+1:j),'r.','MarkerSize',10);
    end
    ylabel('FRET fraction');
    ylim([0,max(ave_FRET)*1.3]);
    %     if max(ave_FRET)>.06
    %         ylim([0,max(ave_FRET)*1.1]);
    %     end
    xlabel(xlabel_input);
    
    drawnow;
    F(j) = getframe(fig1);
    clf;
end
implay(F);
set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 800 600]);
movie2avi(F,'C:\Users\Bryan\Documents\MATLAB\movies\9_9_stage_acclamation','compression','none');
toc