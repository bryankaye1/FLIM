
function [jmax, nph_mean, pixel_counts, dataout] = make_int_bins(ni,num_int_bins,pmat)

%This sub-section combines pixel groups from with similar intensities into SUPER PIXELS
%Number of photons per super pixel is saved in nph_counts. The nph_counts
%field in the input file only exists when this section is run

%Average number of photons per pixel in the super pixel is saved under nph_mean->ni.

intbin = min(ni)+(0:num_int_bins)*(max(ni)-min(ni))/num_int_bins;
nph_mean = [];
dataout = [];
pixel_counts = [];
for k = 1:num_int_bins
    nph = 0;
    datag = 0;
    count =0;
    for l = 1:length(ni)
        if intbin(k) <= ni(l) && intbin(k+1)> ni(l)
            count = count  +1;
            nph = ni(l)+nph;
            datag = pmat(l,:)+datag;
        elseif k==num_int_bins && intbin(k) <= ni(l) && intbin(k+1)>= ni(l)
            count = count  +1;
            nph = ni(l)+nph;
            datag = pmat(l,:)+datag;
        end
    end
    if count > 0
        nph_mean =[nph_mean nph/count];
        dataout = [dataout;datag];
        pixel_counts = [pixel_counts count];     
    end
end
jmax = length(nph_mean);
end