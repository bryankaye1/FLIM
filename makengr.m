function [dataout,ngr,nphi,pixs_per_bin] = makengr(data,pmt,ngr) 
 
[~,hislen] = size(data); %move hislen into spc_2_his_256
%ind is index of pixels, sorted from lowest to highest photon counts
[~, ind] = sort(pmt);

%Here we make dsort, which is a list of histograms, ordered from
%smallest to largest number of photons in that histogram
dsort= pi*ones(length(ind),hislen);
for m = 1:length(ind)
    dsort(m,:) = data(ind(m),:);
    pmtsort(m) = pmt(ind(m));
end
clear data pmt
nphi = pi*ones(ngr,1);
dataout = pi*ones(ngr,hislen);
pixs_per_bin = floor(length(ind)/ngr);
for gr = 1:ngr
    datag = 0;
    pmtgroup = 0;
    %create pixel groups sorted by intensity, with equal number of pixels
    for n = 1 + (gr-1)*pixs_per_bin: gr*pixs_per_bin
        datag = dsort(n,:)+datag;
        pmtgroup = pmtsort(n) + pmtgroup;
    end
    dataout(gr,:) = datag;
    nphi(gr) = pmtgroup;
end

end