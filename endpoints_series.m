%%change nph_counts to pixel_counts after I rerun the int_binned analysis
%In the fit function, I added in a regularization to the weighting to
%(heuristically) account for xvalue uncertainty

%close all;
clear;
uf_sp_start = 12382;%10628;%10568;%
uf_sp_end = 12389;%10687;%10627;%
show_fit=1;
show_images=1;
acq_time = 100;

matstart = uf_sp_start;
matend = uf_sp_end;

%%Colors and marker size for plotting
color1 = [27,158,119]/255; %purple for high acc
color2 = [117,112,179]/255; %green for lower
color3 = [.6,.6,.6]; %Grey for fit lines
%%Load in data, fit model, plot best-fits, save fit parameters
tic
for sample = 1:length(matstart)
    ind = 0;
    for i = matstart(sample):matend(sample)
        clear x y stdpr pixel_counts
        ind = ind + 1;
        %Load in matout files
        [~,output,~] = load_mat_data(i,1);
        imagedata = load_int_image(i);
        %%Set fit paramters for FF v int over image
        ni = output(1,1,1).ni;
        al = output(1,1,1).w1Best./output(1,1,1).w2Best;
        try output(1,1,1).pixel_counts;
            pixc = 1;
            pixel_counts = output(1,1,1).pixel_counts;
        catch
            pixc = 0;
        end
        %%Get FRET fraction and intensity
        for j = 1:length(ni)%length(output(1,1,:))
            x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
                +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
            [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
                output(1,1,j).prestx,al,output(1,1,j).w02est,...
                output(1,1,j).w02estx,x(j),'combine_background');
        end
        %convert from n_photons to cps
        if pixc
            x = x.*(128*128/acq_time); % x = x.*(128*128/20); % x = ph/pixel  .* pixels/sec
        else
            delta_t = 20 / length(ni);
            x = x./delta_t;
        end
        [pfm(sample,ind),Nmon(sample,ind)] = fit_pf_nm(x,y,stdpr,al);
        
        for fill_ind=output(1,1,1).jmax+1:16
            x(fill_ind) = 0;
            y(fill_ind) = 0;
            stdpr(fill_ind) = 1e-5;
            pixel_counts(fill_ind)=0;
        end
        
        xmat(sample,ind,:) = x;
        ymat(sample,ind,:) = y;
        stdmat(sample,ind,:) = stdpr;
        raw_images(:,:,ind) = squeeze(imagedata);
        ti{ind} = output(1,1,1).dataname;
        pf = 0.12;
        
        if pixc
            pol(sample,ind) = sum( (x.*((20/16384).*pixel_counts).*y)...
                ./( pf *( 1 + (al-1).*y) ));
            mon(sample,ind)= sum(x.*((20/16384).*pixel_counts)) - pol(sample,ind)*(1-pf+al*pf);
            frac_pol = pol(sample,ind) / ( pol(sample,ind) + mon(sample,ind) );
            ave_FRET(sample,ind) = sum(y.*pixel_counts./16384);
            nph_bin(ind,:) = pixel_counts.*x.*(acq_time/16384);
        else
            pol(sample,ind) = sum(x.*y./( pf *( 1 + (al-1).*y) ));
            ave_FRET(sample,ind) = mean(y);
            mon(sample,ind)= sum(x) - pol(sample,ind)*(1-pf+al*pf);
            frac_pol(sample,ind) = pol(sample,ind) ./ ( pol(sample,ind) + mon(sample,ind) );
        end
    end
    blackval = min(min(min(raw_images(:,:,2:end))));
    whiteval = max(max(max(raw_images(:,:,2:end))));
    for j = 1:ind
        imseq(:,:,1,j) = mat2gray(raw_images(:,:,j),[blackval whiteval]);
    end
end
%%
xmax = max(max(max(xmat)));
ymax = max(max(max(ymat)));

if show_fit
    indmax = ind;
    F(indmax) = struct('cdata',[],'colormap',[]);
    xmodel = 0:xmax*.001:xmax*1.2;
    fig1 = figure(2); clf;
    for j = 1:indmax
        for sample = 1:length(matstart)
            %subplot(length(matstart),1,sample); hold on;
            subplot(2,2,1); hold on;
            axis([0 xmax*1.2 0 ymax*1.2]);
            errorbar(xmat(sample,j,:),ymat(sample,j,:)...
                ,stdmat(sample,j,:),'.','Color',color1);
            modelb= pf_nm_pred(al,pfm(sample,j),Nmon(sample,j),xmodel);
            %modelb= fitmod(pfm(sample,j),Nmon(sample,j),xmodel);
            plot(xmodel,modelb,'--','Color',color3);
            title(['FF v int fit'],'Interpreter','Tex');
            xlabel('intensity (cps)');
            ylabel('FRET fraction');
            subplot(2,2,[2,4]);
            imshow(imseq(:,:,1,j));
            title(strrep(ti{j},'_',' '));
            subplot(2,2,3); hold on;
            yyaxis left;
            plot(1:length(pol(1:j)),pol(1:j),'.','MarkerSize',10);
            ylabel('eNpol');
            yyaxis right;
            plot(1:length(ave_FRET(1:j)),ave_FRET(1:j),'.','MarkerSize',10);
            ylabel('FRET fraction');
            xlabel('sample #');
            
        end
        drawnow;
        F(j) = getframe(fig1);
        clf;
    end
    implay(F); 
    set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 800 600]);
end
toc
% movie(fig1,F,1,5);
%     if ispc
%         movie_dir = 'C:\Users\Bryan\Desktop\';
%     else
%         movie_dir = '/Users/bryankaye/Desktop';
%     end
%     if output(1,1,1).tfw==0
%         movie_name = [movie_dir,'no_time_removal_6-23-16'];
%     else
%         movie_name = [movie_dir,'\6-23-16'];
%     end
%     movie2avi(F,movie_name,'compression','none');
% end
