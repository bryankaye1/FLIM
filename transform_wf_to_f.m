function [fBest,stdpr, varargout] = transform_wf_to_f(pwf,wfx,al,pwnf,wnfx,x,varargin)
%% Comments
%This function transforms wf (FRET photon fractions) to f(FRET population
% fraction) via bootstrapping
% See supplement for math for transformation
% Threshold for minimum number of photons to try bootstrapping is 50.
%
% Inputs:
% To assume background photons are the same as non-FRET photons:
% varargin{1} = 'combine_background'
% wf = fret photon fraction, wnf = non-fret photon fraction
%pwf: probability of wf aka wf posterior (prest)
%wfx: x-axis of pwf (prestx)
%al: relative brightness of FRETters to non-FRETters, aka alpha. (w1/w2)
%pwnf: probability of wnf aka wnf posterior: (w02est)
%wnfx: x-axis of pwnf (w02estx)
%x: number of photons for bootstrapping (ni)
%
%%
numvarargs = length(varargin);
optargs = {'no_input',40000,'mean','no_input'};
optargs(1:numvarargs) = varargin;
[combine_wnf_wb,N_samples,use_histogram,ig_count] = optargs{:};

if x>500 || strcmp(ig_count,'ignore_pcount')
    if strcmp(combine_wnf_wb,'combine_background')
        wfsams = randsample(wfx,N_samples,true,pwf);
        fsams = wfsams./(wfsams+al*(1-wfsams));
        fsams(isnan(fsams)) = 0;
        stdpr = std(fsams);
    else
        wfsams = randsample(wfx,N_samples,true,pwf);
        wnfsams = randsample(wnfx,N_samples,true,pwnf);
        fsams = wfsams./(wfsams+al*wnfsams);
        fsams(isnan(fsams)) = 0;
        stdpr = std(fsams);
    end
    [yhist,xhist]= histcounts(fsams,200);
    varargout{1} = yhist;
    varargout{2} = xhist;
    if strcmp(use_histogram,'mode')
        [~,max_ind] = max(yhist);
        fBest = xhist(max_ind);
    elseif strcmp(use_histogram,'no_histogram')
        fBest = mode(fsams);
    else
        fBest = mean(fsams);
    end   
else 
    stdpr = 10;
    fBest = 0;
    fprintf('photons less than 500\n');
end

if isnan(stdpr) || isnan(fBest)
    stdpr = 10;
    fBest = 0;
    fprintf('isnan in FF or FF error!\n');
end

%%This section added on 2016-08-10
sample_res = wfx(2)./(wfx(2)+al*(1-wfx(2))) - wfx(1)./(wfx(1)+al*(1-wfx(1)));
if stdpr < sample_res
    fprintf('standard deviation < 2X PDF resolution\n');
end
stdpr = stdpr + sample_res/2;



end

