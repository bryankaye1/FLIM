function [x,y,stdpr] = getFF_int(output,al)
%%Returns FRET fraction and photons  from FRETters.

%If using pixels are divided into SUPER pixels (equal spaced intensities)
%then x is photons per pixel from FRETters

ni = output(1,1,1).ni;
for j = 1:length(ni)
    x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
        +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
    
    if isfield(output,'pixel_counts')
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine',40000,'mean','ignore_pcount');
    else
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine');
    end
end
end