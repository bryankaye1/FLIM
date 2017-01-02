%%Used to plot changes in eNmon with acquisition time

clear;
set(0, 'DefaulttextInterpreter', 'none');
imin = 3384; imax = 3393;

ind = 0;
%%Colors and marker size for plotting
msize = 30;% Marker size
color1 = [27,158,119]/255;%purple for high acc
color2 = [117,112,179]/255;%green for lower
color3 = [.6,.6,.6];%Grey for fit lines
color4 = [0,0,0]; %black

%%Load in data, fit model, plot best-fits, save fit parameters
for i = imin:imax
    %Load in matout files
    clear output intb ffP cintpr
    try
        nstr = strcat('matout\matout',num2str(i),'.mat');
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
            x(j) = (ni(j)/ep(j))*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
               +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
            [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
                output(1,1,j).prestx,al,output(1,1,j).w02est,output(1,1,j).w02estx,x(j));
    end
    %Remove pixel groups with zero photons    
    if min(x)==0
        istart = find(x==0,1,'last')+1;
        x = x(istart:end);
        y = y(istart:end);
        stdpr = stdpr(istart:end);
    end
    
    x = x/(128*128/32);
    %Fit FRET fraction (fluoro pop) v intensity.
    fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
    fresult = fit(x',y',fitmod,'StartPoint',[max(y)/al,min(x)],'Weights',(1./stdpr).^2',...
        'Lower',[0,0],'Upper',[1,max(x)]);
    ab = fresult.a; bmonb = fresult.b;
    ci95 = confint(fresult,.954);
    ci68 = confint(fresult,.682);
    xmodel =0:1000;
    modelb = fitmod(ab,bmonb,xmodel);

    %FF v int plot
    if i ==3385 || i == 3393
        figure(ind); clf; hold on;
        ti = ['Exposure Time: ',num2str(5*ind),' seconds']; title(ti);
        xlabel('Intensity'); ylabel('FRET-fraction');
        h2=errorbar(x,y,stdpr,'.','Color',color2);
        plot(x,y,'.','MarkerSize',msize,'Color',color2);
        plot(xmodel,modelb,'--','Color',color3);
        axis([0 1000 0 .35]);
        set(gca,'XTick',0:200:1000,'YTick',0:.1:.3);
    end
    
    %Save parameters from fits        
    pfm(ind) = ab;%Best Pf
    pferr(ind) = ci95(2,1)-ci95(1,1);% 95% Confidence bound-for errbars
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;% Stdev estimate from confindence int-for fitting weights
    Nmon(ind) = bmonb; %Best Nmon
    nmoner(ind) = ci95(2,2)-ci95(1,2);%95% Confidence bound-for errorbars
    nmonstd(ind) = (ci68(2,2)-ci68(1,2))/2; % Stdev estimate from confindence int-for fitting weights
    %Not Used
    bpol(ind) = sum(x.*y./(ab*al));
    tub(ind) = sum(x);
end

%%
%Plot changes in Nmon
x_data = 5:5:50;
figure(ind+1); clf; hold all;
f = fit(x_data',Nmon','poly1','Weights',(1./nmonstd).^2'); %%Linear fit w/ wieghts
xfit = [0,x_data];
fitvec = f.p1.*xfit+f.p2;
plot(xfit,fitvec,'--','Color',color3);
plot(x_data,Nmon,'.','MarkerSize',msize,'Color',color4);
errorbar(x_data,Nmon,nmoner/2,'.','Color',color3); % This plots 95% condifence bounds
axis([0 55 0 80]);
ax =gca;
ax.XTick = 0:10:50;
ax.YTick= 0:20:85;
fprintf('eNmon vs acq time fit: slope is %f and offset is %f\n',f.p1,f.p2)

%Plot Changes in Nmon
figure(ind+2); clf; hold all;
errorbar(x_data,pfm,pferr/2,'.','Color','k') % This plots 95% condifence bounds
pf_ave = sum(pfm./(pfstd.^2))/sum(1./pfstd.^2);
plot(x_data,pf_ave*ones(1,length(x_data)));
ylabel('Prob FRET');
xlabel('Acquisition Time');
axis([0 55 0 .4]);
ax =gca;
ax.XTick = 0:10:50;
ax.YTick= 0:.1:.4;

