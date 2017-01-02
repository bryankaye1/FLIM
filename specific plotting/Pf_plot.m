%%Used to plot acceptor series plots

clear;
imin = 3384; imax = 3393; %98%84-93


ind = 0;
%%Colors and marker size for plotting
color1 = [27,158,119]/255; %purple for high acc
color2 = [117,112,179]/255; %green for lower
color3 = [.6,.6,.6]; %Grey for fit lines
msize = 30; % Marker size

%%Load in data, fit model, plot best-fits, save fit parameters
for i = 3522:3527    
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
    ep = output(1,1,1).sinti;
    %%Get FRET fraction (fluoro population) and intensity 
    for j = 1:length(output)
            x(j) = (ni(j)/ep(j));%*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
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
    x = x/( 20 * (128*128)/16 );
    %Fit FRET fraction (fluoro pop) v intensity
    fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
    fresult = fit((x)',y',fitmod,'StartPoint',[max(y),min(x)],'Weights',(1./stdpr).^2',...
        'Lower',[0,0],'Upper',[1,max(x)]);
    ab = fresult.a; bmonb = fresult.b;
    ci95 = confint(fresult,.954);
    ci68 = confint(fresult,.682);
    xmodel = 0:400;
    modelb = fitmod(ab,bmonb,xmodel);
    
    %FF v int plots 
    if i == 3524 % Plot of medium acc (green)
       figure(1); clf; hold on
       ax=gca;
       axis([0 400 0 .11]);
       ylim([0 .11])
       ax.XTick = 0:100:400;
       ax.YTick= 0:.02:.1;
       h3=errorbar(x,y,stdpr,'.','Color',color1);
       plot(x,y,'.','MarkerSize',msize,'Color',color1);
       plot(xmodel,modelb,'--','Color',color3);
    elseif i==3526 % Plot of high acc (purple), same figure as green
        figure(1);
        h5=errorbar(x,y,stdpr,'.','Color',color2);
        plot(x,y,'.','MarkerSize',msize,'Color',color2);
        plot(xmodel,modelb,'--','Color',color3);

    elseif i==3527 %Plot Highest acceptor (black), seperate figure 
        figure(2); clf; hold on;
        ax=gca;
        h2=errorbar(x,y,stdpr,'.','Color',color3);
        plot(x,y,'.','MarkerSize',msize,'Color','k');
        plot(xmodel,modelb,'--','Color',color3);
        ax.XTick = 100:100:300;
        ax.YTick= 0:.04:.2; %Was ax.YTick= 0:.03:.12; in 6/2016 submission
        axis([15 350 0 .2]); % Was axis([15 350 0 .12]) in 6/20 submission
        
        %Generate plot to compare to dyes only in extract
        figure(10); clf; hold on;
        ax=gca;
        h2=errorbar(x,y,stdpr,'.','Color',color3);
        plot(x,y,'.','MarkerSize',20,'Color','k');
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
    %Not Used
end

%%

%Plot changes in Pf
figure(3); clf; hold all;
ax = gca();
x = (0:length(pfm)-1);
f = fit(x',pfm','poly1','Weights',(1./pfstd).^2'); %%Linear fit w/ wieghts
fitvec = f.p1.*x+f.p2; %Line with best fit parameters from above
plot(x,fitvec,'--','Color',color3);
plot(x,pfm,'.','MarkerSize',msize,'Color','k');

ylabel('Prob Fret');
xlabel('Sample number');
errorbar(x,pfm,pferr/2,'.','Color','k')
xlim([0 5.5])
ylim([0 .13])
ax.XTick = 0:2.5:5;
ax.YTick= 0:.02:.12;

%Plot Changes in Nmon
figure(4); clf; hold all; ax = gca();
errorbar(x,Nmon,nmoner/2,'.','Color','k') %error
Nmon_ave = sum(Nmon./(nmonstd.^2))/sum(1./nmonstd.^2); %weightd least squares fit - matlab fit fn cannot fit a constant function
plot(x,Nmon_ave*ones(1,length(x)));
ylabel('eNmon');
xlabel('Acceptor Concentration');
axis([-0.25 5.5 0 40]);
ax.XTick = 0:2.5:5;
ax.YTick= 0:10:40;

x = (0:length(pfm)-1)*(1.6/5);
f = fit(x',pfm','poly1','Weights',(1./pfstd).^2'); %%Linear fit w/ wieghts
egg_tub = (18+5+1.6);
nearby_fret_err = confint(f,.954)*egg_tub;
nearby_fret = f.p1*egg_tub;
fprintf('FRETable neighbors are %f +/-%f \n',nearby_fret,nearby_fret-nearby_fret_err(1));
fprintf('Pf vs acc fit: slope is %f and offset is %f\n',f.p1,f.p2);