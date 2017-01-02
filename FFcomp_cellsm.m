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

sept6_2e5 = [17120:17269,17494:17543];
sept9 = 17544:17692;
sept13 = [18092:18141,18042:18091,18142:18191,17992:18041,17841:17891,17942:17991,17892:17941];
sept14_old_irf = [18192:18241,18242:18291,18292:18341];
sept14 = [18342:18391,18492:18541,18392:18441,18542:18591,18442:18491];

%sept21_cells_ngr16 = [18634:18649];

%sept21_w1sweep_ngr16 = [18687,18724,18761,18798,18835,18872...]

w1slice = fix_w1(22244:22280,.4,22354); %fix_w1(first_group,w1_des,last_mat)


ivec = [22555,22100,w1slice];%[22133:22243,21726:21762]; %[jul31_donor_ran,jul31_uf1_ran,jul31_uf2_ran];
acq_time = 10;
%%%
show_fit=1;
show_images=1;
save_time = acq_time+10;
start_time = 120;
total_pixels = 128*128;
num_pho_thresh=5000;
tic
ind = 0;
for i = ivec
    clear x y stdpr pixel_counts
    ind = ind + 1;
    %Load in data and images
    [~,output,flagoutput] = load_mat_data(i,1);
    if flagoutput
        pausehere=1;
    end
    
    imagedata = load_int_image(round(i));
    %%Get FRET fraction and intensity. For now we make x(j)=ni
    ni = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    for j = 1:length(ni);
        x(j) = ni(j);%*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
        %+sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_background');
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
        if isfield(output,'pixel_counts')
            [pfm(ind),Nmon(ind)] = fit_pf_nm(x_int,y,stdpr,al);
        else
            [pfm(ind),Nmon(ind)] = fit_pf_nm(x,y,stdpr,al);
        end
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
    tbac = num2str(output(1,1,1).tbac);
    tfw = num2str(output(1,1,1).tfw);
    w1 = num2str(output(1,1,1).w1min);
    w2 = num2str(output(1,1,1).w2min);
    ti{ind} = [strrep(strrep(output(1,1,1).dataname,'.sdt',''),'_',' '),...
        ' trem:',tbac,'/',tfw, ' w1/w2: ', w1,'/',w2];
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
    x_plot = (start_time:save_time:start_time+save_time*(ind-1))./60;
    x_limits = [x_plot(1)-save_time/2, x_plot(end)+save_time/2];
    xlabel_input = 'minutes';
    
    
    
end

sumxmat = sum(xmat,2);
F(ind) = struct('cdata',[],'colormap',[]);
fig1 = figure(2); clf;
for j = 1:ind
    %Place fit in figure (or total photons
    
    if show_fit && length(ni)>1
        fig1= figure(3);
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
        axis([0 xmax*1.2 0 .1]);
        title(ti{j});
    else %Or just total number of photons
        % plot(1:j,sumxmat(1:j)','go');
        % ylabel('Total Photons'); %xlim(x_limits);
    end
    
    %Place intensity image in figure
    %     pos_im = [.65 .3 .4 .4];
    %     subplot('Position',pos_im);
    %     imshow(imseq(:,:,1,j));
    %     if isempty(c1index)
    %     title(strrep(ti{j},'_',' '));
    %     else
    %     base_name = strrep(ti{j},'_',' ');
    %     ctindex = strfind(strrep(ti{j},'0',''),'_c');
    %     base_name = base_name(1:ctindex);
    %     tmin = num2str(floor(x_plot(j)));
    %     tsec = num2str(round(mod(x_plot(j),1)*60),'%02.0f');
    %     title([base_name,' ',tmin,':',tsec]);
    %     xlabel([num2str(ivec(1)),'-',num2str(ivec(end))])
    %     end
    %
    %     %Plot total photons, eNpol and mean FRET fraction
    %     if show_fit && length(ni)>1
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
    F(j) = getframe(fig1);
    clf;
end
implay(F);
set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 800 600]);
%movie2avi(F,'/Users/bryankaye/Desktop/thaw','compression','none');
toc