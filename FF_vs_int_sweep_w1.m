%%Used to plot acceptor series plots

clear;
matstart = 8011;%8581;
matend = 8067;%8637;%

%%Colors and marker size for plotting
color1 = [27,158,119]/255; %purple for high acc
color2 = [117,112,179]/255; %green for lower
color3 = [.6,.6,.6]; %Grey for fit lines
msize = 10; % Marker size
set(0, 'DefaulttextInterpreter', 'none')
ind = 0;
%%Load in data, fit model, plot best-fits, save fit parameters


for i = matstart:matend
    ind = ind + 1;
    %Load in matout files
    clear output intb ffP cintpr x y stdpr
    
    try
        nstr = strcat('Y:\Users\bkaye\cluster\matout\matout',num2str(i),'.mat');
        load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
    catch
        fprintf('\nmatout %s DNE',num2str(i));
        break;
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
    xmodel(ind,:) = 0:max(x)*.01:max(x)*1.2;
    modelb(ind,:) = fitmod(ab,bmonb,xmodel(ind,:));
    xmovie(ind,:) = x;
    ymovie(ind,:) = y;
    errmovie(ind,:) = stdpr;
    w1movie(ind) = output(1,1,1).w1max;
    
    %Save parameters from fits
    pfm(ind) = ab; %Best Pf
    pferr(ind) = ci95(2,1)-ci95(1,1);% 95% Confidence bound-for errbars
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;% Stdev estimate from confindence int-for fitting weights
    Nmon(ind) = bmonb; %Best Nmon
    nmoner(ind) = ci95(2,2)-ci95(1,2);%95% Confidence bound-for errorbars
    nmonstd(ind) = (ci68(2,2)-ci68(1,2))/2;% Stdev estimate from confindence int-for fitting weights
    pol(ind) = sum(x.*y./( ab *( 1 + (al-1).*y) ));
    ave_FRET(ind) = mean(y);
    
end
%%
indmax = ind;
F(indmax) = struct('cdata',[],'colormap',[]);
%cps = 1e6;
%delta_t = ( sum(sum(xmovie)) /indmax ) / cps;
delta_t = 100 / length(ni);
xmax = max(max(xmovie))/delta_t;
ymax = max(max(ymovie));
for j = 1:indmax
    
    figure(1); hold on;
    axis([0 xmax*1.2 0 ymax*1.2]);
    errorbar(xmovie(j,:)./delta_t,ymovie(j,:),errmovie(j,:),'.','Color',color1);
    plot(xmodel(j,:)./delta_t,modelb(j,:),'--','Color',color3);
    title([strrep(output(1,1,1).dataname,'_',' '),' - \tau_{short} is ',...
        num2str(w1movie(j)),'(ns)'],'Interpreter','Tex');
    xlabel('intensity (cps)');
    ylabel('FRET fraction');
    drawnow
    F(j) = getframe(gcf);
    clf;
end
fig = figure(2);
movie(fig,F,1,5);
if output(1,1,1).tfw==0
movie_name = [strrep(output(1,1,1).dataname,'_',' '),' no time removal'];
else
movie_name = [strrep(output(1,1,1).dataname,'_',' '),' time removal'];  
end
movie2avi(F,movie_name,'compression','none');
%%
figure(1); clf; 
xvec = .2:.05:3;
yyaxis left; hold on;
ylabel('Pf and mean FRET fraction');
plot(xvec,pfm);
plot(xvec,ave_FRET);
yyaxis right;
plot(xvec,pol);
ylabel('eNpol');
xlabel('\tau_{short} (ns)','Interpreter','Tex');
legend('Pf','mean FRET level','eNpol','Location','North');
title(strrep(output(1,1,1).dataname,'_',' '),'Interpreter','Tex');

%%
