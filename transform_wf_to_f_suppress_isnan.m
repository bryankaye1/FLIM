function [fBest,stdpr,varargout] = transform_wf_to_f_suppress_isnan(pwf,wfx,al,pwnf,wnfx,x)

if x>50
    wfsams = randsample(wfx,40000,true,pwf);
    wnfsams = randsample(wnfx,40000,true,pwnf);
    fsams = wfsams./(wfsams+al*wnfsams);
    fsams(isnan(fsams)) = 0;
    fBest = mean(fsams);
    stdpr = std(fsams);
     
    [yhist,xhist]= hist(fsams,50);
    
    varargout{1} = yhist;
    varargout{2} = xhist;
    
else
    stdpr = 10;
    fBest = 0;
    fprintf('photons less than 50\n');
end

if isnan(stdpr) || isnan(fBest)
    stdpr = 10;
    fBest = 0;
    fprintf('isnan in FF or FF error!\n');
end

end

