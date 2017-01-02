% Enter "type" = 1 to read out posterior of lifetime
% Enter "type" = 0 to read matouts to plot # of short lifetimes photons vs longlifetime photons
% enter "int" = low/med/hi for auto plots

clear;
imin = 2744; imax = 2744;
type=9;
ind =0;
for i=imin:imax
    
    try
        nstr = strcat('matout',num2str(i),'.mat');
        tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
    catch exception
        try
            nstr = strcat('Z:\bkaye\cluster\mof\matout',num2str(i),'.mat');
            tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
        catch exception
            continue
        end
    end
    ind =ind+1;
    output = tempf.output(:,:,:);% output = output(20,2,2);
    jmax = tempf.output(1,1,1).jmax;
    exptmax = tempf.output(1,1,1).exptmax;
    cyclesmax = tempf.output(1,1,1).cyclesmax;
    
    for jind = 1:jmax
        for expt = 1:exptmax
            for cindex = 1:cyclesmax
                
                w2Best(cindex,expt,jind) = output(cindex,expt,jind).w2Best;
                w1Best(cindex,expt,jind) = output(cindex,expt,jind).w1Best;
                prBest(cindex,expt,jind) = output(cindex,expt,jind).prBest;
                w02Best(cindex,expt,jind) = output(cindex,expt,jind).w02Best;
                
                %  cintpr(cindex,expt,jind) = findci(output(cindex,expt,jind).prest,...
                % output(cindex,expt,jind).prestx,.682);
                
                %  figure(2*cindex); clf; plot(output(cindex,expt,jind).w1estxl1,output(cindex,expt,jind).w1estl1);
                %  figure(2*cindex-1); clf; plot(output(cindex,expt,jind).w1estx,output(cindex,expt,jind).w1est);
                %  figure(2*cindex); clf; plot(output(cindex,expt,jind).prestxl1,output(cindex,expt,jind).prestl1);
                %  figure(2*cindex-1); clf; plot(output(cindex,expt,jind).prestx,output(cindex,expt,jind).prest);
            end
        end
    end
    prB(ind)= mean(prBest);
    if type==3
        poli = prBest*w2Best/w1Best;
        pol(mod(ind-1,6)+1,floor((ind-1)/6)+1) = poli;
        fprintf('%1.1f,%1.1f is %s\n', mod(ind-1,6)+1, floor((ind-1)/6)+1,output(1,1,1).dataname);
    end
    
    
    % prBtemp(ind) = mean(prBest);
    if type==0 || type==2
        
        w1B = mean(w1Best);
        w2B = mean(w2Best);
        w02B = mean(w02Best);
        
        prB(ind)= mean(prBest);
        w01(ind,:) = prBest.*w02B.*w1B./(w2B.*(1-prBest));
        w01B(ind) = mean(w01(ind,:));
        w01std(ind) = std(w01(ind,:));
        int = output(1,1,1).dataname(end-5:end);
        
    else if type==1
            prB(ind)= mean(prBest);
            figure(ind); clf; plot(output(1,1,1).w1estx,output(1,1,1).w1est);
            fprintf('%f for %s\n', w1Best, output(1,1,1).dataname);
            w2B(ind) = w1Best;
            if ind == imax-imin+1
                w2m = mean(w2B);
                w2std = std(w2B);
                fprintf('mean/std is %1.4f/%1.4f\n', w2m,w2std);
            end
        end
        
        
    end
end
%%


if type==4
    figure(13); clf; hold all;
    %plot(length(pfm),pfm,'-o');
    f = fit((0:length(prB)-1)',prB','poly1');
    plot(f,0:length(prB)-1,prB,'o');
    axis([0 length(prB) 0 max(prB)*1.2]);
    ylabel('FRET Frac');
    xlabel('Sample number');
end


if type==2
    
    fest = w01B;
    if strcmp(output(1,1,1).dataname(11:14),1024)
        iratio = str2double(output(1,1,1).dataname(11:14));
    else
        iratio = str2double(output(1,1,1).dataname(11:13));
    end
    ratio = 1/iratio;
    fact = ratio*ones(1,length(fest));
    w01se = w01std;%/sqrt(cyclesmax);
    
    figure(1); clf; hold all;
    plot(fest,'bo');
    plot(fact,'r');
    
    errorbar(1:length(fest),fest,w01se,w01se,':','Color',[.7 .7 .9]); %errorbar(X,Y,L,U)
    ti1 = sprintf('Est of short fract vs Pho Number \n %s control dye short frac = 1/%s',int, num2str(iratio,3));
    title(ti1);
    xlabel('Number of Photons in Millions'); % x-axis label
    ylabel('Standard Deviation'); % y-axis label
    axis([0,9,0,max(fest+w01std)])
    set(gca,'XTickLabel',['    '; '50.0';'25.0';'12.5';'6.25';'3.13';'1.56';'0.78';'0.39']);
    
    
    figure(2); clf; hold all;
    plot(log2(w01std),'bo');
    plot(.5*(0:7)+log2(w01std(1)),'--');
    ti2 = sprintf('Std Dev vs Pho Number \n %s control dye short frac = 1/%s',int, num2str(iratio,3));
    set(gca,'XTickLabel',['    '; '50.0';'25.0';'12.5';'6.25';'3.13';'1.56';'0.78';'0.39']);
    title(ti2);
    xlabel('NPho is 50M/2^x'); % x-axis label
    ylabel('Standard Deviation'); % y-axis label
    
end
%%
if type ==0
    ratio = 1./2.^(1:ind);
    fest = w01B;
    fact = ratio;
    wlu = log2((w01std/sqrt(cyclesmax) + fest)./fest);
    wll = log2(fest./(fest-w01std/sqrt(cyclesmax)));
    
    if strcmp(int,'low')
        figure(2);
    elseif strcmp(int,'med')
        figure(4);
    elseif strcmp(int,'hi')
        figure(6);
    end
    
    clf; hold all;
    plot(log2(fest),'bo');
    plot(log2(fact),'r');
    
    errorbar(1:length(fest),log2(fest),wll,wlu,':','Color',[.7 .7 .9]); %errorbar(X,Y,L,U)
    ti = strcat(int,' int - fret frac ');%strcat('tfw/tbac = ',num2str(output(1,1,1).tfw),' / ', num2str(output(1,1,1).tbac));%, '- w2 is ', num2str(w2B));
    title(ti);
    
    clear wll wlu
    for j = 1:ind-1
        
        s1 = w01(j,:);
        s2 = w01(j+1,:);
        s4 =[];
        for k = 1:cyclesmax-1
            s3 = s1 - circshift(s2,k,2);
            s4 = [s4 s3];
        end
        fechse(j) = std(s4)/sqrt(cyclesmax);
        fech(j) = fest(j)-fest(j+1);
        if fech(j)<0
            fech(j) =0;
            fechse(j) =0;
        end
    end
    
    
    for m = 1:ind-1
        if fech(m)-fechse(m) < 0
            wll(m) = 4;
            wlu(m) = log2((fechse(m)+fech(m))/fech(m));
            fprintf('spot %f: errorbar goes below zero\n', m');
        elseif fech(m)==0
            wll(m) = 0;
            wlu(m) = 0;
            fprintf('spot %f: fret change positive\n', m');
        else
            wlu(m) = log2((fechse(m)+fech(m))/fech(m));
            wll(m) = log2(fech(m)/(fech(m)-fechse(m)));
        end
    end
    
    if strcmp(int,'low')
        figure(1);
    elseif strcmp(int,'med')
        figure(3);
    elseif strcmp(int,'hi')
        figure(5);
    end
    
    clf; hold all;
    plot(log2(fech),'bo');
    plot(log2(fact(2:end)),'r');
    errorbar(1:length(fech),log2(fech),wll,wlu,':','Color',[.7 .7 .9]);
    ti2=strcat(int,' int - Changes in fret frac');
    title(ti2);
end
%%
if type ==3
    figure(1); clf;
    mpol = mean(pol,1);
    stdpol = std(pol,0,1)/sqrt(2);
    polsub = mpol(2:end) - mpol(1);
    errorbar(0:4,mpol,stdpol(1:end));
end
%%

% if type ==0
% l = 5*10^7;
% ratio = 1./2.^(1:18);
% s = l*ratio;
% w1 = .4572;
% w2 = 3.93;
%
% prpred = s .*(w2/w1) ./ (s.*(w2/w1)+l);
% figure(2); clf;   hold all;
% plot(log(prB),'b');
% plot(log(prpred),'r');
%
% wlu = log((prstd + prB)./prB);
% wll = log((prstd - prB)./prB);
%
% errorbar(1:length(prB),log(prB),wll,wlu,':','Color',[.7 .7 .9]); %errorbar(X,Y,L,U)
% ti = strcat('tfw/tbac = ',num2str(output(1,1,1).tfw),' / ', num2str(output(1,1,1).tbac));%, '- w2 is ', num2str(w2B));
% title(ti);
% end





% mean(w1Best)
% std(w1Best)
% mean(prBest)
% std(prBest)
% tfw = 0:.1:1;

%plot(tfw,mw1);
% hold all;
% errorbar(tfw,mw1,sw1,'-o');
% %xlabel(

%plot(tfw,mpr);

% rpr = [.1,.1,.1,.1,.1,.1];
% rw1 = [.5,1,1.5,.5,1,1.5];
%
% figure(1);hold all;
% scatter([rpr(1),rpr(4)],[mpr(1),mpr(4)]);
% scatter([rpr(2),rpr(5)],[mpr(2),mpr(5)]);
% scatter([rpr(3),rpr(6)],[mpr(3),mpr(6)]);
%
% errorbar(.1,.1,max(sep));
% title('Est w1 vs actual (simulated) w1');
% legend('w1=.5','1ns','1.5ns','Stdev','Location','NorthEast');
%
% figure(2); hold all; scatter(rw1,mw1); %errorbar(rw1,mw1,sw1);
% title('Est FF vs actual (simulated) FF');
% errorbar(rw1(1:3),rw1(1:3),sw1(1:3));


%%
% figure(1); plot(log2(rpr),log2(rpr),'--gs',...
%     'LineWidth',0.1,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);
% %
% % figure(2); plot(rpr,rpr,'+');
%
% %%
% figure(2);
% E = mean(spr,1);
% errorbar(rpr,rpr,E);

% rpr = [.1,.1,.1,.1,.1,.1];
% rw1 = [.5,1,1.5,.5,1,1.5];
%
% figure(1); plot(rpr,mpr);
% figure(2); plot(rpr,mpr);


% figure(1); errorbar(log2(rpr),log2(rpr),log2(E));

% figure(3); clf; hold all;
% y1 =  squeeze(mean(w2Bestmat(:,1,:),1));
% y2 =  squeeze(mean(w2Bestmat(:,2,:),1));
% x = .1:.1:1;
% plot(x,y1); plot(x,y2);

% for i = 1:10
% y1(:) = squeeze(w2Bestmat(:,1,i))
% y2(:) = squeeze(w2Bestmat(:,2,i))
% x = [i/10, i/10, i/10, i/10, i/10];
% scatter(x,y1);
% scatter(x,y2);
%
% end

% title('2014-6-13 Donor only + Ran. pr set to 0');
% xlabel('time removed after peak (ns)');
% ylabel('Non-FRET (w2) Lifetime');
% legend('irf 235pm', 'irf 642pm','Location','SouthEast');


%%from comment section in for loop
%     a = [a; w1Best];
%     dataname = output.dataname;
%     mw1 = mean(squeeze(w1Best));
%     sw1 = std(squeeze(w1Best));
%
%     mw2 = mean(squeeze(w2Best));
%     sw2 = std(squeeze(w2Best));
%     fprintf('w1 is %f\n',w1Best(6));
%     fprintf('DN is %s\n',output(6,1).dataname);

% outc = struct2cell(output);
%errorv = squeeze(cell2mat(outc(60,:,:,:)));

%     w1Best(ind) = squeeze(cell2mat(outc(59,:,:,:)));
%     w2Best = squeeze(cell2mat(outc(58,:,:,:)));
%     prBest(ind)= squeeze(cell2mat(outc(57,:,:,:)));

%     mpr(ind) = mean(prBest);
%     spr(ind) = std(prBest);
%     sep(ind) = std(prBest)/sqrt(length(prBest));
%
%     mw1(ind) = mean(w1Best);
%     sw1(ind) = std(w1Best);
%     sw1(ind) = std(w1Best)/sqrt(length(w1Best));
%
%comment = output(1,1,1).comment; %fprintf('%s\n\n',comment);

% filenames = outc(26,:,:,:);
% if sum(errorv)==0
%fprintf('There are no errors\n');
% else
%     fprintf('There is an error in %f\n',ind);
% end


% fprintf('prmean is %f prpred is %f \n', mean(prBest), prpred(ind));


%     pol = 4*prBest;
%     cint = 4*cintpr;
%     ti = output(1,1).dataname;
%     figure(ind); clf;% scatter(1:length(prBest),prBest,'w');
%     hold all;
%    title(ti(1:end-8));% title('Time-Series 2uM Filtered Extract'); %
%     ylabel('Estimated Polymer fraction');
%     xlabel('Time Point');
%     %axis([0,30,.06,.2]);
%     errorbar(pol,cint,':','Color',[.7 .7 .9]);
%     plot(1:length(pol),pol,'r');
%     legend('68% Conf int','Post Max', 'Location','NorthEast');


