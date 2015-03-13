function [data,jmax, nphi,varargout] = spc_2_his(tmini,tmaxi,file_name,pth_sdt,ngr,tsm,varargin)
%sort pixels by intensity, equal pixels per group


for ts=1:tsm
    if length(file_name)>3 %fixes filename if sdt is appended to file name
        if strcmp(file_name(end-3:end),'.sdt')
            file_name=file_name(1:end-3);
        end
    end
    if tsm ==1
        file_name2 = strcat(file_name,'.sdt');
    else
        file_name2 = strcat(file_name,'_c',num2str(ts),'.sdt');
    end
    sdt = bh_readsetup([pth_sdt file_name2]); block=1;
    ch = bh_getdatablock(sdt,block); %raw lifetime data
    ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
    ld = double(ld);  
    if length(size(ld))==3
        datat(:,:,1+128*(ts-1):128*ts) = ld(tmini:tmaxi,:,:);
        hislen = tmaxi-tmini+1;
    end
end

if nargin > 6
    sdt = bh_readsetup([varargin{1} varargin{2}]);     block=1; 
    ch = bh_getdatablock(sdt,block); %raw lifetime data
    imap1 = double(ch); %raw lifetime data in single form (normally stored as 4 bit number)   
    imap1 = imap1(tmini:tmaxi,:,:);
    imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image   
else
    imap1 = pi*ones(4096,128,128); %%If you disnt input int map, we use a pi map
    imap1 = imap1(tmini:tmaxi,:,:);
    imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image    
end    
    if  length(size(ld))==2 %checks to see if this is an image (3dimensional) or just a FIFO (2 dimensional)
        dataout = ld(tmini:tmaxi);
        nphi = sum(dataout);
        ngr = 1;
    elseif ngr==1
        dataout = squeeze(sum(sum(datat,3),2));
        nphi = sum(dataout);
    else
        ivec = pi*ones(128*128*tsm,1);%
        pmt = pi*ones(128*128*tsm,1);%
        datat2 = pi*ones(hislen,128*128*tsm);
        for k = 1:128*tsm
            for j=1:128           
                ivec(j+(k-1)*128) = sum(imap(:,j,k)); %total photon number, reshaped into line
                pmt(j+(k-1)*128) = sum(datat(:,j,k)); %total photon number, reshaped into line
                datat2(:,j+(k-1)*128) = datat(:,j,k); %flimage, linearized
            end
        end
        [~, ind] = sort(pmt); %i is index of pixels, low to high pho counts
       
        %Here we make dsort, which is a list of histograms, ordered from
        %smallest to largest number of photons in that histogram
        dsort= pi*ones(hislen,length(ind));
        for m = 1:length(ind)
            dsort(:,m) = datat2(:,ind(m));
        end
        
        %%Here we ensure that the number of groups is such that every
        %%photon belongs to one group
        while mod(length(ind),ngr)~=0
            ngr = ngr+1;
        end
        
        sinti = pi*ones(1,ngr);
        nphi = pi*ones(1,ngr);
        dataout = pi*ones(hislen,ngr);
        for gr = 1:ngr
            datag = 0;
            sint = 0;
            for n = (1+(gr-1)*length(ind)/ngr):1:gr*length(ind)/ngr
                datag = dsort(:,n)+datag;
                sint = ivec(n)+sint; %sum of ints
            end
            dataout(:,gr) = datag;
            nphi(gr) = sum(datag);
            sinti(gr) = sint;
        end
    end
    varargout{1} = sinti;
    data=dataout';
    jmax = ngr; %g-1; % REPLACE "g-1" in other thrt = 2,3 in code  
end