function [x,y,stdpr] = getFF_int(output,al,varargin)
%%Returns FRET fraction and photons  from FRETters.

%If using pixels are divided into SUPER pixels (equal spaced intensities)
%then x is photons per pixel from FRETters

if length(varargin)
    n_samples = varargin{1};
else
    n_samples = 40000;
end

ni = output(1,1,1).ni;
for j = 1:length(ni)
    x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
        +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
    
    if isfield(output,'pixel_counts')
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine',n_samples,'mean','ignore_pcount');
    else
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine',n_samples);
    end
end
end