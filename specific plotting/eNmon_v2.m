clear;
imin = 3384; imax = 3393; %98%84-93
figure(1); clf; hold on; 
color1 = [.6,.6,.6];
color2 = [27,158,119]/255;%[217,95,2]/255;
color3 = [117,112,179]/255;
msize = 30;

ind = 0;
for i = 3384:3393%[3385,3393]   %[3384,3389,3393]
    clear output intb ffP cintpr
    try
        nstr = strcat('matout',num2str(i),'.mat');
        tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
    catch exception
        try
            nstr = strcat('Y:\Users\bkaye\cluster\mof\matout',num2str(i),'.mat');
            tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
        catch
            continue
        end
    end
    ind =ind+1;
    output = tempf.output(:,:,:);
    dataname = output.dataname;
    ni = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    ep = output(1,1,1).sinti;
    for j = 1:length(output)
        
        x(j) = (ni(j)/ep(j))*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
            +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,output(1,1,j).w02estx,x(j));
    end
    
    if min(x)==0
        istart = find(x==0,1,'last')+1;
        x = x(istart:end);
        y = y(istart:end);
        stdpr = stdpr(istart:end);
    end
    
    fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
    fresult = fit(x',y',fitmod,'StartPoint',[max(y)/al,min(x)],'Weights',(1./stdpr).^2',...
        'Lower',[0,0],'Upper',[1,max(x)]);
    ab = fresult.a; bmonb(ind) = fresult.b; 
    ci = confint(fresult,.68);
    nmoner(ind) = ci(2,2)-ci(1,2);
    
    xmmax = 5*10^5;
    xm = 100:100:xmmax;
    modelb = fitmod(ab,bmonb(ind),xm);
    
    if i==3385
        axis([0 xmmax 0 .35]);
        h1=errorbar(x,y,stdpr,'.','MarkerSize',msize,'Color',color2);
        plot(xm(modelb>0),modelb(modelb>0),'--','Color',color1);
    elseif i== 3393
        h2=errorbar(x,y,stdpr,'.','MarkerSize',msize,'Color',color3);
        plot(xm(modelb>0),modelb(modelb>0),'--k','Color',color1);

    end

end

ax =gca;
ax.XTick = (1:1:5)*10^5;
ax.YTick= 0:.1:.35;
%ax.XTickLabel = {'','','',''};
%ax.YTickLabel={'','','',''};

% legend([h1,h2,h3],{'5sec:  Pf= 34%','30sec: Pf= 32%','50sec: Pf= 33%'},'Location','SouthEast');


%%
plotype='Nmon';
if strcmp(plotype,'Nmon')
    x = 5:5:50;
    f = fit(x',bmonb','poly1');
    xfit = [0,x];
    fitvec = f.p1.*xfit+f.p2;
    
    figure(ind+1); clf; hold all;
    ax =gca;
    %ax.XTick = (1:1:5)*10^5;
    ax.YTick= 0:10^4:4*10^4;
    errorbar(x,bmonb,nmoner,'.','Color','k');
    plot(xfit,fitvec,'--','Color',color1);
    plot(x,bmonb,'.','MarkerSize',msize,'Color','k');
    axis([0 max(x)*1.1 0 max(bmonb)*1.1]);
    xlabel('Acquisition Time (s)');
    ylabel('\epsilon Nmon');
end

% old =0;
% if old ==1
%
%     bpolr = reshape(bpol,3,5);
%     mbpol = mean(bpolr,1);
%     bpolmin = max(mbpol)/2^(length(mbpol)-1);
%     intx = bpolmin*2.^(0:length(mbpol)-1);
%     figure(ind+1); clf; hold all;
%     plot(intx,mbpol,'bo'); plot(intx,intx,'r'); %errorbar(intx,intm,intv);
%     title('Pol by FRET');
%
%     tubr = reshape(tub,3,5);
%     mtub = mean(tubr,1);
%     tubmin = max(mtub)/2^(length(mtub)-1);
%     tubx = tubmin*2.^(0:length(mtub)-1);
%
%     figure(ind+2); clf; hold all;
%     plot(tubx,mtub,'bo'); plot(tubx,tubx,'r'); %errorbar(intx,intm,intv);
%     title('polymer by intensity');
%
%     figure(ind+3); clf; hold all; plot(tubx,mbpol,'bo'); plot(tubx,tubx,'r');
%     title('polymer by FRET vs Polymer by intensity');
% end
