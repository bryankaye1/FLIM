%%change nph_counts to pixel_counts after I rerun the int_binned analysis
%In the fit function, I added in a regularization to the weighting to
%(heuristically) account for xvalue uncertainty

%close all;
clear;
matstart = 10688;%10628;%10568;%
matend = 10747;%10687;%10627;%
acq_time = 20;
%%%
show_fit=1;
show_images=1;
save_time = acq_time+5;
start_time = 45;
total_pixels = 128*128;
num_pho_thresh=5000;
tic
ind = 0;
for i = matstart:matend
    clear x y stdpr pixel_counts
    ind = ind + 1;
    %Load in matout files

    [~,output,~] = load_mat_data(i,1);
    imagedata = load_int_image(round(i));

    %%Set fit paramters for FF v int over image
    ni = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    try pixel_counts = output(1,1,1).pixel_counts;
        pixc = 1;
    catch
        pixc = 0;
    end
    %%Get FRET fraction and intensity
    for j = 1:length(ni);
        x(j) = ni(j);%*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
            %+sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'combine_background');
    end
    %convert from n_photons to cps
    if pixc
        x_int = x.*(total_pixels/acq_time); % x = x.*(128*128/20); % x = ph/pixel  .* pixels/sec
        x = pixel_counts.*x_int.*(acq_time/total_pixels);
    else
        delta_t = acq_time / length(ni);
        x_int = x./delta_t;
    end
    if show_fit && length(ni)>1
        [pfm(ind),Nmon(ind)] = fit_pf_nm(x,y,stdpr,al);
    end
    if pixc
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
    ti{ind} = strrep(output(1,1,1).dataname,'.sdt','');
    pf = 0.12;
    
    if pixc
    pol(ind) = sum( x.*y./( pf *( 1 + (al-1).*y) ));
    ave_FRET(ind) = sum(y.*x./(1-y+al*y))/sum(x./(1-y+al*y));     
        
        %  nph_bin(ind,:) = pixel_counts.*x.*(acq_time/total_pixels);
       % pol(ind) = sum( nph_bin(ind,:).*y./( pf *( 1 + (al-1).*y) ));
       % mon(ind)= sum(nph_bin(ind,:)) - pol(ind)*(1-pf+al*pf);
       % frac_pol = pol(ind) / ( pol(ind) + mon(ind) );
      %  ave_FRET(ind) = sum(y.*pixel_counts./total_pixels);     
    else
        pol(ind) = sum(x.*y./( pf *( 1 + (al-1).*y) ));
       % mon(ind)= sum(x) - pol(ind)*(1-pf+al*pf);
       % frac_pol(ind) = pol(ind) ./ ( pol(ind) + mon(ind) );
        ave_FRET(ind) = mean(y);
    end
    
    thr_pho{ind} = find(xmat(ind,:)<num_pho_thresh);
end

blackval = min(min(min(raw_images(:,:,1:end))));
whiteval = max(max(max(raw_images(:,:,1:end))));
for j = 1:ind
    imseq(:,:,1,j) = mat2gray(raw_images(:,:,j),[blackval whiteval]);
end
%%
c1index = strfind(strrep(ti{1},'0',''),'_c');
if isempty(c1index)
    x_plot = 1:ind;
    x_limits = [.5,ind+.5];
    xlabel_input = 'sample number';
else
    x_plot = (start_time:save_time:45+save_time*(ind-1))./60;
    x_limits = [x_plot(1)-save_time/2, x_plot(end)+save_time/2];
    xlabel_input = 'minutes';
end

if pixc
    sumnbin = sum(xmat,2);
   % sumnbin = sum(nph_bin,2);
else
    sumxmat = sum(xmat,2);
end

F(ind) = struct('cdata',[],'colormap',[]);
fig1 = figure(2); clf;
for j = 1:ind
    %Place fit in figure (or total photons
    subplot(2,2,1);
    if show_fit && length(ni)>1
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
    else
        if pixc  
        plot(1:j,sumnbin(1:j)','o');    
        else
        plot(1:j,sumxmat(1:j)','o');
        end
        ylabel('Total Photons'); xlim(x_limits);
    end
    
    %Place intensity image in figure
    pos_im = [.5 .3 .5 .5];
    subplot('Position',pos_im);
    imshow(imseq(:,:,1,j));
    if isempty(c1index)
    title(strrep(ti{j},'_',' '));
    else
    base_name = strrep(ti{j},'_',' ');
    base_name = base_name(1:c1index);
    tmin = num2str(floor(x_plot(j)));
    tsec = num2str(round(mod(x_plot(j),1)*60),'%02.0f');   
    title([base_name,' ',tmin,':',tsec]);    
    xlabel([num2str(matstart),'-',num2str(matend)])
    end
    
    subplot(2,2,3); hold on;
    yyaxis left;
    plot(x_plot(1:j),pol(1:j),'.','MarkerSize',10);
    if pixc
        plot(x_plot(1:j),sumnbin(1:j)'*pol(1)/sumnbin(1));
    else
        plot(x_plot(1:j),sumxmat(1:j)'/pol(1)/sumxmat(1));
    end
    ylabel('eNpol & Total Photons');
    yyaxis right;
    plot(x_plot(1:j),ave_FRET(1:j),'.','MarkerSize',10);
    ylabel('FRET fraction');
    xlabel(xlabel_input);
    
    drawnow;
    F(j) = getframe(fig1);
    clf;
end
implay(F);
set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 800 600]);
toc