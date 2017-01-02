%%Used to plot acceptor series plots

clear;
imin = 3384; imax = 3393; %98%84-93


ind = 0;
%%Colors and marker size for plotting
color1 = [27,158,119]/255; %purple for high acc
color2 = [117,112,179]/255; %green for lower
color3 = [.6,.6,.6]; %Grey for fit lines
msize = 10; % Marker size

%%Load in data, fit model, plot best-fits, save fit parameters
for i = [26212,25955]
    %Load in matout files
    clear output intb ffP cintpr
    try
        nstr = strcat('matout',num2str(i),'.mat');
        tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
    catch exception
        try
            nstr = strcat('Y:\Users\bkaye\cluster\matout\matout',num2str(i),'.mat');
            tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
        catch
            continue
        end
    end
    ind =ind+1;
    output = tempf.output(:,:,:);
    dataname = output.dataname;
    %%Set fit paramters for FF v int over image
    ni = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    %ep = output(1,1,1).sinti;
    %%Get FRET fraction (fluoro population) and intensity
    for j = 1:length(output)
        x(j) = ni(j);%(ni(j)/ep(j));%*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
        %  +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,output(1,1,j).w02estx,x(j));
    end
    %Remove pixel groups with zero photons
    if min(x)==0
        istart = find(x==0,1,'last')+1;
        x = x(istart:end);
        y = y(istart:end);
        cintpr = cintpr(istart:end);
    end
    int_norm = 86/33;
    x = ( x/( 20 * (128*128)/16 ) ) * int_norm;
    %Fit FRET fraction (fluoro pop) v intensity
    fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
    fresult = fit((x)',y',fitmod,'StartPoint',[max(y),min(x)],'Weights',(1./stdpr).^2',...
        'Lower',[0,0],'Upper',[1,max(x)]);
    ab = fresult.a; bmonb = fresult.b;
    ci95 = confint(fresult,.954);
    ci68 = confint(fresult,.682);
    xmodel = 0:400;
    modelb = fitmod(ab,bmonb,xmodel);
    
    if i == 26212
    figure; clf; hold on;
    ax=gca;
    h2=errorbar(x,y,stdpr,'.','Color',color3);
    plot(x,y,'.','MarkerSize',msize,'Color','k');
    plot(xmodel,modelb,'--','Color',color3);
    ax.XTick = 100:100:300;
    ax.YTick= 0:.03:.12; %Was ax.YTick= 0:.03:.12; in 6/2016 submission
    axis([15 350 0 .12]); % Was axis([15 350 0 .12]) in 6/20 submission
    end
    
    %Save parameters from fits
    pfm(ind) = ab; %Best Pf
    pferr(ind) = ci95(2,1)-ci95(1,1);% 95% Confidence bound-for errbars
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;% Stdev estimate from confindence int-for fitting weights
    Nmon(ind) = bmonb; %Best Nmon
    nmoner(ind) = ci95(2,2)-ci95(1,2);%95% Confidence bound-for errorbars
    nmonstd(ind) = (ci68(2,2)-ci68(1,2))/2;% Stdev estimate from confindence int-for fitting weights
end
