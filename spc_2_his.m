function [data,jmax, nphi,varargout] = spc_2_his(tmini,tmaxi,file_name,pth_sdt,ngr,tsm,varargin)
%sort pixels by intensity, equal pixels per group
%set tsm for 1 exposure. tsm > 1 for time series
% nphi is number of photons per group
%jmax is how many bins there are
%data is flim decay curve
%adjusted photon count per pixel group



%FOV intensity correction,
%Vargargin(1) = FOV intensity map file path
%Vargargin(2) = FOV intensity map file name

%Intensity image, add 'imout; to last varargin
%Varargin(1) = 'imout' for image out and with no FOV correction
%Varargin(3) = 'imout' for image out with FOV intensity correction

%varagout{1} = image out if varargin1 or varargin3 is 'imout'


%ts/tsm is for combining multiple exposures into one image
%concontanates flimages into "datat"

%changes 3/21/16: changed code to correct for intensity at the (not grouped)
%pixel level.
for ts=1:tsm
    if length(file_name)>3 %fixes filename if sdt is appended to file name
        if strcmp(file_name(end-3:end),'.sdt')
            file_name=file_name(1:end-4);
        end
    end
    if tsm ==1
        file_name2 = strcat(file_name,'.sdt');
    else
        if tsm <10
            file_name2 = strcat(file_name,'_c',num2str(ts),'.sdt');
        elseif tsm < 100
            if ts<10
                file_name2 = strcat(file_name,'_c0',num2str(ts),'.sdt');
            else
                file_name2 = strcat(file_name,'_c',num2str(ts),'.sdt');
            end
        end
        
    end
    ld = sdt_to_vector(pth_sdt,file_name2);
    if length(size(ld))==3
        datat(:,:,1+128*(ts-1):128*ts) = ld(tmini:tmaxi,:,:);
        hislen = tmaxi-tmini+1;
    end
end

%%In this section we seperate the data into different groups, and keep
%%tract of intensity

%First check if varargins wer passed
if nargin > 6
    %If imout was passed as first varargout, return an intensity image of
    %the data (last frame for multiple exposures).
    if strcmp(varargin{1},'imout')
        varargout{1} = sum(ld,1);
        % If file path for the intensity map was passed (via varargout 1 and
        % 2), create intensity map
        
    elseif strcmp(varargin{1},'FLIMage')
        if nargin > 8
            pth_int = varargin{2};
            fn_int = varargin{3};
            if length(fn_int)>3 %fixes filename if sdt is appended to file name
                if strcmp(fn_int(end-3:end),'.sdt')
                    fn_int=fn_int(1:end-4);
                end
            end
            fn_int = strcat(fn_int,'.sdt');
            imap1 = sdt_to_vector(pth_int, fn_int);
            imap1 = imap1(tmini:tmaxi,:,:);
            imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image
        else
            ivec = 1;
            imap1 = pi*ones(4096,128,128); %%If you do not input int map, we use a pi map
            imap1 = imap1(tmini:tmaxi,:,:);
            imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image
        end

        for k = 1:128
            for j = 1:128               
                ivec(j+(k-1)*128) = sum(imap(:,j,k)); %total photon number of intensity map, reshaped into line
                pmt(j+(k-1)*128) = sum(datat(:,j,k)); %total photon number of data, reshaped into line                
                dataout(:,j+(k-1)*128)=datat(:,j,k);
            end
        end       
        data = dataout';
        jmax = 128*128;
        nphi = pmt./ivec;
        varargout{1} = pmt;
        return
    else
        pth_int = varargin{1};
        fn_int = varargin{2};
        if length(fn_int)>3 %fixes filename if sdt is appended to file name
            if strcmp(fn_int(end-3:end),'.sdt')
                fn_int=fn_int(1:end-4);
            end
        end
        fn_int = strcat(fn_int,'.sdt');
        imap1 = sdt_to_vector(pth_int, fn_int);
        imap1 = imap1(tmini:tmaxi,:,:);
        imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image
        
        %If imout was passed as third varargout, return an intensity image of
        %the data (last frame for multiple exposures).
%         if exist(varargin{3},'var') 
%             if strcmp(varargin{3},'imout')
%             varargout{2} = sum(ld,1);
%             end
%         end
    end
else
    imap1 = pi*ones(4096,128,128); %%If you do not input int map, we use a pi map
    imap1 = imap1(tmini:tmaxi,:,:);
    imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image
end

%% return FLIM vectors
%If not an image, return a single vector
if  length(size(ld))==2
    dataout = ld(tmini:tmaxi);
    nphi = sum(dataout);
    ngr = 1;
    varargout{1} = 1;
    %if ngr is 1, group pixels from image into FLIM vector
elseif ngr==1
    dataout = squeeze(sum(sum(datat,3),2));
    nphi = sum(dataout);
    varargout{1} = 1;
else
    %Otherwise build grouped pixels, number of photons in each pixel group (ni),
    %and relative intensity (sinti)
    ivec = pi*ones(128*128*tsm,1);%
    pmt = pi*ones(128*128*tsm,1);%
    datat2 = pi*ones(hislen,128*128*tsm);
    for k = 1:128*tsm
        for j=1:128
            ivec(j+(k-1)*128) = sum(imap(:,j,k)); %total photon number of intensity map, reshaped into line
            pmt(j+(k-1)*128) = sum(datat(:,j,k)); %total photon number of data, reshaped into line
            datat2(:,j+(k-1)*128) = datat(:,j,k); %flimage, linearized
        end
    end
    
    %Here we renormalize pixel intensity before sorting
    pmt = pmt./(ivec/mean(ivec));    
    %ind is index of pixels, sorted from lowest to highest photon counts
    [~, ind] = sort(pmt);
    
    %Here we make dsort, which is a list of histograms, ordered from
    %smallest to largest number of photons in that histogram
    dsort= pi*ones(hislen,length(ind));
    for m = 1:length(ind)
        dsort(:,m) = datat2(:,ind(m));
        %int_map_sort(m) = ivec(ind(m));
        pmtsort(m) = pmt(ind(m));
    end
    %int_map_sort = int_map_sort/mean(int_map_sort);
    
    %%Here we ensure that the number of groups is such that every
    %%photon belongs to one group
    while mod(length(ind),ngr)~=0
        ngr = ngr+1;
    end
    
%    sinti = pi*ones(1,ngr);
    nphi = pi*ones(1,ngr);
    dataout = pi*ones(hislen,ngr);
    for gr = 1:ngr
        datag = 0;
%        sint = 0;
        pmtgroup = 0;
        for n = (1+(gr-1)*length(ind)/ngr):1:gr*length(ind)/ngr
            datag = dsort(:,n)+datag;
            %sint = pmtsort(n)/int_map_sort(n)+sint; %sum of ints
            pmtgroup = pmtsort(n) + pmtgroup;
        end
        dataout(:,gr) = datag;
        nphi(gr) = sum(pmtgroup);
        %sinti(gr) = pmtgroup/sint; %This devision returns effective int_map_sort (or epsilon)
    end
    %varargout{1} = sinti;
end
data=dataout';
jmax = ngr; %g-1; % REPLACE "g-1" in other thrt = 2,3 in code
end