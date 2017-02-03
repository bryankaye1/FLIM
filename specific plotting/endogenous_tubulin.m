%%change nph_counts to pixel_counts after I rerun the int_binned analysis
%In the fit function, I added in a regularization to the weighting to
%(heuristically) account for xvalue uncertainty

%close all;
clear;
color3 = [.6,.6,.6]; %Grey for fit lines
msize = 30; %size of marker
jan10 = [27226,27228:27231,27250,27232,27233,27227,27251,27234:27237];
jan10_paperlife = [27252,27254:27259,27253,27265,27260:27264];
ivec = jan10;%27252:27265;%[27252,27254:27257];

%%%
acq_time = 400;

ind = 0;
dataname_ind = {};
for i = ivec
    clear x y stdpr pixel_counts
    ind = ind + 1;
    %Load in data and images
    [~,output,flagoutput] = load_mat_data(i,1);

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
    
    delta_t = acq_time / length(ni);
    x_int = x./delta_t;
    x = x_int;
    [pfm(ind),Nmon(ind),fresult] = fit_pf_nm(x,y,stdpr,al);
    ci68 = confint(fresult,.682);
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;
    
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
    thr_pho{ind} = find(xmat(ind,:)<10000);
    ti{ind} = strrep(output(1,1,1).dataname,'.sdt','');
    
end
% for j = 1:ind
%     figure(j); clf;
%     xmax = max(max(xintmat));
%     ymax = max(max(ymat));
%     hold on;
%     errorbar(xintmat(j,:),ymat(j,:),stdmat(j,:),'.');
%     errorbar(xintmat(j,thr_pho{j}),ymat(j,thr_pho{j}),stdmat(j,thr_pho{j}),'.','Color','red');
%     xmodel = 0:xmax*.001:xmax*1.2;
%     plot(xmodel,pf_nm_pred(al,pfm(j),Nmon(j),xmodel),'--','Color',[0.6,0.6,0.6]);
%     title(ti{j},'Interpreter','Tex');
%     xlabel('intensity (cps)');
%     ylabel('FRET fraction');
%     axis([0 xmax*1.2 0 ymax*1.2]);
% end



%%
figure(ind+1); clf; hold all;
% ax = gca();
% x = (0:length(pfm(1:4))-1);
% f = fit(x',pfm(1:4)','poly1','Weights',(1./pfstd(1:4)).^2'); %%Linear fit w/ wieghts
% fitvec = f.p1.*x+f.p2; %Line with best fit parameters from above
% plot(x,fitvec,'--','Color',color3);
% plot(x,pfm(1:4),'.','MarkerSize',msize,'Color','k');

x = (0:length(pfm)-1);
f = fit(x',pfm','poly1','Weights',(1./pfstd).^2'); %%Linear fit w/ wieghts
fitvec = f.p1.*x+f.p2; %Line with best fit parameters from above
plot(x,fitvec,'--','Color',color3);
plot(x,pfm,'.','MarkerSize',msize,'Color','k');

xlabel('sample number');
ylabel('Pf');

%%
m = 3*0.074638*10;
b = f.p2;


a = 2*16.8/(140+2+1.75);
d = 1.75*25/(140+2+1.75);
ul_stock = 30;
u_vol = 0:4;
u = ul_stock.*u_vol./(21+u_vol);
Pf = pfm(end-4:end);

err = pfstd(end-4:end);
x=u;
y=Pf;



fitmod = @(m,e,x)(m*a./(a+d+x+e) +b );
fit_e = fit((x)',y',fitmod,'StartPoint',[0.074638/10,1.6],'Weights',...
    (1./(err)).^2','Lower',[0,0],'Upper',[100,50]);

figure(2);
plot(x,y,'.','MarkerSize',msize,'Color','k');
hold on; 
xmodel = 0:.01:max(x);
plot(xmodel,fitmod(fit_e.m,fit_e.e,xmodel),'--','Color',color3);


%e_vec = (-a*b -b*d -a*m +a*Pf +d*Pf -b*u +Pf.*u)./(b -Pf);



%%

%
%
% %Determine if this is a time series and set appropriate plot parameters
% c1index = strfind(strrep(ti{1},'0',''),'_c');
% if isempty(c1index)
%     x_plot = 1:ind;
%     x_limits = [.5,ind+.5];
%     xlabel_input = 'sample number';
%     subpvec = []
% else
%     [start_nums,end_nums,~] = find_series(dataname_ind);
%     x_plot = [];
%     %spacer = [60*5,1000,60*8,0]./tpoint_dur;
%     spacer = 15;
%     for i1 = 1:length(start_nums)
%     x_plot = [x_plot,start_time:tpoint_dur:...
%         start_time+tpoint_dur*(end_nums(i1)-start_nums(i1))];
%     start_time = start_time+tpoint_dur*(end_nums(i1)-start_nums(i1)+spacer+1);
%     end
%     x_plot = x_plot./60;
%     x_limits = [x_plot(1)-tpoint_dur/2, x_plot(end)+tpoint_dur/2];
%     xlabel_input = 'minutes';
% end
%
% sumxmat = sum(xmat,2);
% F(ind) = struct('cdata',[],'colormap',[]);
% fig1 = figure(2);
% for j = 1:ind
%     clf;
%     %Place fit in figure (or total photons
%     if show_fit && length(ni)>1
%         subplot(2,2,1);
%         xmax = max(max(xintmat));
%         ymax = max(max(ymat));
%         hold on;
%         errorbar(xintmat(j,:),ymat(j,:),stdmat(j,:),'.');
%         errorbar(xintmat(j,thr_pho{j}),ymat(j,thr_pho{j}),stdmat(j,thr_pho{j}),'.','Color','red');
%         xmodel = 0:xmax*.001:xmax*1.2;
%         plot(xmodel,pf_nm_pred(al,pfm(j),Nmon(j),xmodel),'--','Color',[0.6,0.6,0.6]);
%         title('FF v int fit','Interpreter','Tex');
%         xlabel('intensity (cps)');
%         ylabel('FRET fraction');
%         axis([0 xmax*1.2 0 ymax*1.2]);
%     else %Or just total number of photons
%        % plot(1:j,sumxmat(1:j)','go');
%        % ylabel('Total Photons'); %xlim(x_limits);
%     end
%
%     %Place intensity image in figure
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
%
%     drawnow;
%     F(j) = getframe(fig1);
%
% end
% implay(F);
% set(findall(0,'tag','spcui_scope_framework'),'position',[50 50 1200 1000]);
% %movie2avi(F,'C:\Users\Bryan\Documents\MATLAB\movies\oct29_30minc_noran','compression','none');
% toc

