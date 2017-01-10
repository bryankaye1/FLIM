 function [model,varargout] = fitFLIM(output,qv,ti)
%     varargout{1} = pn;
%     varargout{2} = res;
%     varargout{3} = sumres;

    if ispc
    addpath('Y:\Users\bkaye\cluster\matlab_scripts\');
    else
    %addpath(genpath(pwd));
    addpath('./cluster_scripts/matlab_scripts/');
    end 
    
    T=12.58;
    bins = output(qv).bins;
    brem = output(qv).brem;
    ga = output(qv).ga;
    binskeep = bins - brem;
  
    s = T/bins:T/bins:T; %time vector used to generate PDF of signal exponential
    wig = output(qv).wig; %wig = wig';
    tfw =output(qv).tfw; %forward amount of time to remove
    tbac =output(qv).tbac; %backwards amount of time to remove
    p = output(qv).datahis;
    w2 = output(qv).w2Best;
    w02 = output(qv).w02Best;
    w01 = output(qv).prBest;
    [srem,erem,p2,wig2] = remove_bins(T,bins,tfw,tbac,p,wig); %returns data, wigs, and indeces of which bins to keep/remove
    
    f2 = exp(-s/w2); %signal over one period
    f2 = [f2 f2]; %signal over 2 consecutive periods
    f2con = conv(f2,ga); %PDF after conv
    f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
    f2h = f2bar(1:binskeep); %Keep only the appropriate bins
    f2h = f2h/sum(f2h);
    
    back = (1-w01-w02)/bins;
    model = (f2h + back).*wig;
    
    model = model/sum(model);
    pn = p/sum(p);
    res = pn-model;
    sumres = sum(res.^2);
    
    varargout{1} = pn;
    varargout{2} = res;
    varargout{3} = sumres;
     
     figure; clf; hold all;
     plot(log(pn), '.','MarkerSize',5, 'Color', 'b'); plot(log(model),'r'); title(ti);
    
end
    