%%Used to plot acceptor series plots

clear;
matstart = [8011,8068,8125,8182];%,8239,8296,8353,8410];
matend = [8067,8124,8181,8238];%,8295,8352,8409,8466];

%%Colors and marker size for plotting
color1 = [27,158,119]/255; %purple for high acc
color2 = [117,112,179]/255; %green for lower
color3 = [.6,.6,.6]; %Grey for fit lines
msize = 10; % Marker size
%set(0, 'DefaulttextInterpreter', 'none')

%%Load in data, fit model, plot best-fits, save fit parameters
tic
for sample = 1:length(matstart)
    ind = 0;
    for i = matstart(sample):matend(sample)
        ind = ind + 1;
        %Load in matout files
        clear output intb ffP cintpr x y stdpr
        
        try
            nstr = strcat('C:\Users\Bryan\Documents\MATLAB\data\matout\',num2str(i),'.mat');
            load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
        catch
            try
                nstr = strcat('Y:\Users\bkaye\cluster\matout\matout',num2str(i),'.mat');
                load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
            catch
                fprintf('\nmatout %s DNE',num2str(i));
                break;
            end
        end
        
        %%Set fit paramters for FF v int over image
        ni = output(1,1,1).ni;
        
        al = output(1,1,1).w1Best./output(1,1,1).w2Best;
        %ep = output(1,1,1).sinti;
        ex = output(1,1,1).cyclesmax;
        %%Get FRET fraction and intensity
        
        for j = 1:length(ni)%length(output(1,1,:))
            x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
                +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
            
            [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
                output(1,1,j).prestx,al,output(1,1,j).w02est,...
                output(1,1,j).w02estx,x(j),'combine_background');
        end
        delta_t = 100 / length(ni);
        x = x./delta_t;
        %Remove pixel groups with zero photons
        if min(x)==0
            istart = find(x==0,1,'last')+1;
            x = x(istart:end);
            y = y(istart:end);
            cintpr = cintpr(istart:end);
        end
        %Fit FRET fraction (fluoro pop) v intensity
        fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
        fresult = fit((x)',y',fitmod,'StartPoint',[max(y),min(x)],'Weights',(1./stdpr).^2',...
            'Lower',[0,0],'Upper',[1,max(x)]);
        ab = fresult.a; bmonb = fresult.b;
        ci95 = confint(fresult,.954);
        ci68 = confint(fresult,.682);
        
        xmovie(sample,ind,:) = x;
        ymovie(sample,ind,:) = y;
        errmovie(sample,ind,:) = stdpr;
        w1movie(sample,ind) = output(1,1,1).w1max;
        
        %Save parameters from fits
        pfm(sample,ind) = ab; %Best Pf
        pferr(sample,ind) = ci95(2,1)-ci95(1,1);% 95% Confidence bound-for errbars
        pfstd(sample,ind) = (ci68(2,1)-ci68(1,1))/2;% Stdev estimate from confindence int-for fitting weights
        Nmon(sample,ind) = bmonb; %Best Nmon
        nmoner(sample,ind) = ci95(2,2)-ci95(1,2);%95% Confidence bound-for errorbars
        nmonstd(sample,ind) = (ci68(2,2)-ci68(1,2))/2;% Stdev estimate from confindence int-for fitting weights
        pol(sample,ind) = sum(x.*y./( ab *( 1 + (al-1).*y) ));
        ave_FRET(sample,ind) = mean(y);
        
    end
    %ni(sample) = ni;
    dataname{sample} = strrep(output(1,1,1).dataname,'_',' ');
end
toc
%%
indmax = ind;
F(indmax) = struct('cdata',[],'colormap',[]);
%cps = 1e6;
%delta_t = ( sum(sum(xmovie)) /indmax ) / cps;
%delta_t = 100 / length(ni);
xmax = max(max(max(xmovie)));
ymax = max(max(max(ymovie)));
xmodel = 0:xmax*.001:xmax*1.2;
%f_width = 400;
%f_height = 800;

fig1 = figure(1); clf; %'position',[10 10 f_width f_height]); clf;
for j = 1:indmax
    %title(['\tau_{short} is ',num2str(w1movie(j)),'(ns)'],'Interpreter','Tex');
    %axis([0 xmax*1.2 0 ymax*1.2]);
    for sample = 1:length(matstart)
        subplot(4,round(length(matstart)/4),sample); hold on;
        axis([0 xmax*1.2 0 ymax*1.2]);
        errorbar(xmovie(sample,j,:),ymovie(sample,j,:)...
            ,errmovie(sample,j,:),'.','Color',color1);
        modelb= fitmod(pfm(sample,j),Nmon(sample,j),xmodel);
        plot(xmodel,modelb,'--','Color',color3);
        title([dataname{sample},' - \tau_{short} is ',num2str(w1movie(sample,j)),'(ns)'],'Interpreter','Tex');
        xlabel('intensity (cps)');
        ylabel('FRET fraction');
    end
    drawnow;
    F(j) = getframe(fig1);
    clf;
end
%[h, w, p] = size(F(1).cdata);
%fig2 = figure(2);
%set(fig2,'position',[150 150 f_width f_height])
movie(fig1,F,1,5);
if output(1,1,1).tfw==0
    movie_name = 'C:\Users\Bryan\Desktop\16_pixel_groups_no_time_removal_6-23-16';
else
    movie_name = 'C:\Users\Bryan\Desktop\16_pixel_groups_no_time_removal_6-23-16';
end
movie2avi(F,movie_name,'compression','none');
%%
figure(3); clf;
xvec = .2:.05:3;
y1max = max(max(pfm));
y2max = max(max(pol));


for sample = 1:length(matstart)  
    subplot(4,round(length(matstart)/4),sample);
    title(dataname{sample},'Interpreter','Tex');
    
    yyaxis left; hold on;
    ylabel('Pf and mean FRET fraction');
    ylim([0 y1max]);
    plot(xvec,pfm(sample,:));
    plot(xvec,ave_FRET(sample,:));
    
    yyaxis right;
    ylim([0 y2max])
    plot(xvec,pol(sample,:));
    ylabel('eNpol');
    xlabel('\tau_{short} (ns)','Interpreter','Tex');

    legend('Pf','mean FRET level','eNpol','Location','North');
end
%%
