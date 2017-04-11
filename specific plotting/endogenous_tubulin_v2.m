%%change nph_counts to pixel_counts after I rerun the int_binned analysis
%In the fit function, I added in a regularization to the weighting to
%(heuristically) account for xvalue uncertainty

%close all;
clear;
color3 = [.6,.6,.6]; %Grey for fit lines
msize = 30; %size of marker
jan10 = [27226,27228:27231,27250,27232,27233,27227,27251,27234:27237];
jan10_paperlife = [27252,27254:27259,27253,27265,27260:27264];
mar7 = [28339,28340,28341,28344,28342,28343,28346,28345,28347];%[28339,28341,28344,28342,28345,28347]; %
mar16 = [28395:28397,28400,28398,28399,28401:28403];%[28395:28397,28400,28399,28401:28403];%
ivec = mar16;%27252:27265;%[27252,27254:27257];

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
    %al = .45;
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
    [pfm(ind),nmon(ind),fresult] = fit_pf_nm(x,y,stdpr,al);
    ci68 = confint(fresult,.682);
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;
    
    figure(ind); clf; hold on; errorbar(x,y,stdpr,'ko');
    xmodel = 0:max(x)*.001:max(x)*1.2;
    plot(xmodel,pf_nm_pred(al,pfm(ind),nmon(ind),xmodel),'--','Color',[0.6,0.6,0.6]);
    ti{ind} = strrep(output(1,1,1).dataname,'.sdt',''); title(ti{ind});
    
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
    
    if ind ==4
        %bryan =1;
    end
end

%%
a = 2.5*16.8/23;
a35 = a*3/5;
d = 450*2.5/230;

ul_stock = 100;
u_vol = [0,0.5,0.75,1.0,1.5,3];
%u_vol = [0,0.75,1.0,1.5,3];
u = ul_stock.*u_vol./23;
e0 = 16; %between 23.5 and 23.75 is the bifurcation point

num_lin = 3;
yl = pfm(1:num_lin);
stdline = pfstd(1:num_lin);


err = pfstd(num_lin+1:end);
x=u;
y=pfm(num_lin+1:end);

for i = 1:1
    xl = [0,a35/(a35+e0+d),a/(a+e0+d)];
%fit linear regime
    fitline = @(m,bl,x)(m * x + bl);
    f = fit(xl',yl',fitline,'Weights',(1./stdline).^2',...
        'StartPoint',[2,.01]); %%Linear fit w/ wieghts

    %fit unlabeled titration regime
    m =  2;%f.m; %2; f.m*(18+5+1.6);
    b = 0;%f.bl; %0;%

    fitmod = @(e,x)(m*a./(a+x+d+e) + b );
    fit_e = fit((x)',y',fitmod,'StartPoint',[20],'Weights',...
        (1./(err)).^2','Lower',[0],'Upper',[500]);    
%      fitmod = @(m,e,x)(m*a./(a+x+d+e) + b );
%      fit_e = fit((x)',y',fitmod,'StartPoint',[2,20],'Weights',...
%          (1./(err)).^2','Lower',[0,0],'Upper',[10,50]);

    e0 = fit_e.e;
    fprintf('m,e,b are %2.2f,%2.2f,%2.2f\n',m,e0,b);
end

figure(ind+1); clf; hold on;
fitvec = f.m.*xl+f.bl; %Line with best fit parameters from above
plot(xl,fitvec,'--','Color',color3);
errorbar(xl,yl,pfstd(1:num_lin),'.','MarkerSize',msize,'Color','k');
xlabel('% tubulin labeled');
ylabel('Pf');

figure(ind+2); clf;
errorbar(x,y,pfstd(num_lin+1:end),'bo');
hold on;
xmodel = -1:.01:max(x)*1.1;
%plot(xmodel,fitmod(fit_e.m,fit_e.e,xmodel),'--','Color',color3);
plot(xmodel,fitmod(fit_e.e,xmodel),'--','Color',color3);
xlabel('micromolar of added unlabeled tubulin');
ylabel('Pf');
xlim([-1 max(xmodel)]);


figure(12);
y=nmon(num_lin+1:end);
tub = a+u+d+e0;
x = 1./tub;
plot(x,y,'bo');
xlabel('inverse micromolar of added unlabeled tubulin');
ylabel('eNmon');

figure(13); clf; hold on;
y=nmon(num_lin+1:end);
tub = a+u+d+e0;
x = 1./tub;
sam_ind = [1,3,5,6];
x = x(sam_ind);
y = y(sam_ind);
plot(x,y,'bo');
f_nmon = fit(x',y','poly1');
xmodel = linspace(min(x)*.9,max(x)*1.1,1000);
plot(xmodel,f_nmon.p1*xmodel+f_nmon.p2);
xlabel('inverse micromolar of added unlabeled tubulin');
ylabel('eNmon');
%ADD 2D FIT BELOW
%fitobject = fit([x,y],z,fitType)
%%
%ADD 2D FIT HERE
%fitobject = fit([x,y],z,fitType)














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

