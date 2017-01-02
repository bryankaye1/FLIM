%clear; %close all;
imin = 3923; imax = imin; %98%84-93
fprintf('matname is matin%3.0f',imin)
w2=0; w1 =1;

ffbin =0;
isamp = 3;
ind = 0;

fit =1;
for i = imin:1:imax
    clear output intb ffP cintpr
    try
        nstr = strcat('matin',num2str(i),'.mat');
        load(nstr,'-mat','input');
    catch exception
        try
            nstr = strcat('Y:\Users\bkaye\cluster\matin\matin',num2str(i),'.mat');
            load(nstr,'-mat','input');
        catch
            continue
        end
    end
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
        ind = ind+1;
    end
    output = tempf.output;
    w2best = tempf.w2Bestmat;
    
    if w2 ==1
        [m,n,p] = size(output);
        for q = 1:n
            fitFLIM(output,[1,q,1],2*q,sprintf('\nspot-%1.0f w2: %1.3f',q,w2best(q)));
            figure(2*q-1);clf; 
            
            plot(output(1,q,1).w2estx,output(1,q,1).w2est);
            
            title(sprintf('\nspot-%1.0f w2: %1.3f',q,w2best(q)));
            fprintf('\nw2 for spot-%1.0f is: %1.2f',q,w2best(q));
            [ebinary,estring] = read_matout_errorcheck(output(1,q,1).error,output(1,q,1).errparam);
        end
        fprintf('\n\n The mean is : %1.3f\n',mean(w2best));
    end
    
    if w1 ==1
        figure(i-imin+1);
        plot(output(1,1,1).w1estx,output(1,1,1).w1est);
        w1best(i-imin+1) = output(1,1,1).w1Best;
        fprintf('\nw1 for spot-%1.0f is: %1.3f',i-imin+1,w1best(i-imin+1));
        [ebinary,estring] = read_matout_errorcheck(output(1,1,1).error,output(1,1,1).errparam);
    end
    if ffbin ==1
        ffB(ind) = output(1,1,1).prBest; fprintf('\nffB is %1.3f', ffB(ind));
        ff = output(1,1,1).prest;
        ffx = output(1,1,1).prestx;
        figure(i-imin+1);
        plot(ffx,ff);
        [ebinary,estring] = read_matout_errorcheck(output(1,1,1).error,output(1,1,1).errparam);
    end
    
end

if ffbin==1
    ffB = reshape(ffB,isamp,ind/isamp);
    %plot(0:length(mean(ffB,1))-1,mean(ffB,1),'o'); xlabel('Sample #');
    %ylabel('Fret Fraction');
end

if w1 ==1
    fprintf('\n\n w1 mean is: %1.3f\n',mean(w1best));
end
%%
% if fit ==1
%     q=3;
%     qv = [1,q,1];
%     T=12.58;
%     bins = output(qv).bins;
%     brem = output(qv).brem;
%     ga = output(1,q,1).ga;
%     binskeep = bins - brem;
%   
%     s = T/bins:T/bins:T; %time vector used to generate PDF of signal exponential
%     wig = output(qv).wig; wig = wig';
%     tfw =output(qv).tfw; %forward amount of time to remove
%     tbac =output(qv).tbac; %backwards amount of time to remove
%     p = output(qv).datahis;
%     w2 = output(qv).w2Best;
%     w02 = output(qv).w02Best;
%     w01 = output(qv).prBest;
%     [srem,erem,p2,wig2] = remove_bins(T,bins,tfw,tbac,p,wig); %returns data, wigs, and indeces of which bins to keep/remove
%     
%     f2 = exp(-s/w2); %signal over one period
%     f2 = [f2 f2]; %signal over 2 consecutive periods
%     f2con = conv(f2,ga); %PDF after conv
%     f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
%     f2h = f2bar(1:binskeep); %Keep only the appropriate bins
%     f2h = f2h/sum(f2h);
%     
%     back = (1-w01-w02)/bins;
%     model = (f2h + back).*wig;
%     
%     model = model/sum(model);
%     pn = p/sum(p);
%     figure(3); clf; hold all; plot(pn, 'b'); plot(model,'r'); 
%     
%     %res = (data-model)./sqrt(data);
%     %sumresi = sum(res.^2);
%     
%     
%     
%     
% end




