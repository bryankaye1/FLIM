function [data,jmax, nphi,varargout] = spc_2_his_256(tmini,tmaxi,file_name,pth_sdt,ngr,tsm,varargin)
%sort pixels by intensity, equal pixels per group
%set tsm for 1 exposure. tsm > 1 for time series
% nphi is number of photons per group
%jmax is how many bins there are
%data is flim decay curve
%adjusted photon count per pixel group



%FOV intensity correction,
%Varargin(1) = FOV intensity map file path
%Varargin(2) = FOV intensity map file name

%If you want FLIMage that is intensity-map adjusted (for super pixels etc)
%Varargin(1) = 'FLIMout'
%Varargin(2) = FOV intensity map file path
%Varargin(3) = FOV intensity map file name

%Intensity image, add 'imout; to last varargin
%Varargin(1) = 'imout' for image out and with no FOV correction
%Varargin(3) = 'imout' for image out with FOV intensity correction

%varagout{1} = image out if varargin1 or varargin3 is 'imout'


%ts/tsm is for combining multiple exposures into one image
%concontanates flimages into "datat"

%changes 3/21/16: changed code to correct for intensity at the (not grouped)
%pixel level.
if tsm ==0
    tsm = 1;
end
for ts=1:tsm
    file_name = remove_sdt(file_name); %removes ".sdt" from end of filename if present
    if tsm ==1
        file_name2 = strcat(file_name,'.sdt');
    else
        if tsm <10
            file_name2 = strcat(file_name,'_c0',num2str(ts),'.sdt');
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
        ld = make128(ld);
        datat(:,:,1+128*(ts-1):128*ts) = ld(tmini:tmaxi,:,:);
        hislen = tmaxi-tmini+1;
    end
end

%%In this section we seperate the data into different groups, and keep
%%track of intensity

%First check if varargins wer passed
if nargin > 6
    %If imout was passed as first varargout, return an intensity image of
    %the data (last frame for multiple exposures).
    if strcmp(varargin{1},'imout')
        varargout{1} = sum(ld,1);
    elseif strcmp(varargin{1},'save_int_image')
        data = sum(ld,1);
        jmax = 1;
        nphi = pi;
        if nargin > 8
            pth_int = varargin{2};
            fn_int = varargin{3};
            fn_int = remove_sdt(fn_int); %removes ".sdt" from end of filename if present
            fn_int = [fn_int,'.sdt'];
            imap = make128(sdt_to_vector(pth_int, fn_int));
            FOV_int = sum(imap(tmini:tmaxi,:,:),1);
            FOV_int = FOV_int/mean(mean(FOV_int));
            data = squeeze(data)./squeeze(FOV_int);
        end
        return
        
    elseif strcmp(varargin{1},'FLIMage')
        if nargin > 7
            %This section sets the intensity-normalization map
            pth_int = varargin{2};
            fn_int = varargin{3};
            fn_int = remove_sdt(fn_int); %removes ".sdt" from end of filename if present
            fn_int = [fn_int,'.sdt'];
            imap1 = make128(sdt_to_vector(pth_int, fn_int));
            imap1 = imap1(tmini:tmaxi,:,:);
            imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image
        else
            %This section sets the intensity-normalization map to null (if
            %no info is given)
            ivec = 1;
            imap1 = pi*ones(4096,128,128); %%If you do not input int map, we use a pi map
            imap1 = imap1(tmini:tmaxi,:,:);
            imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image
        end

        if nargin > 9
            %This section will boxcar sum pixels together
            reach = varargin{4};
            if reach>0
            datat = boxcar_averager(datat,reach);
            imap = boxcar_averager(imap,reach);
            end
        end
        
        for k = 1:128*tsm
            for j = 1:128               
                ivec(j+(k-1)*128) = sum(imap(:,j,k)); %total photon number of intensity map, reshaped into line
                pmt(j+(k-1)*128) = sum(datat(:,j,k)); %total photon number of data, reshaped into line                
                dataout(:,j+(k-1)*128)=datat(:,j,k);
            end
        end 
        data = dataout';
        jmax = 128*128;
        nphi = pmt./(ivec./mean(ivec));
        return % This makes it so we dont overwrite our output later in the "return FLIM vectors" section
    else
        pth_int = varargin{1};
        fn_int = varargin{2};
        fn_int = remove_sdt(fn_int); %removes ".sdt" from end of filename if present
        fn_int = [fn_int,'.sdt'];
        imap1 = sdt_to_vector(pth_int, fn_int);
        imap1 = make128(imap1);
        imap1 = imap1(tmini:tmaxi,:,:);
        imap = repmat(imap1,1,1,tsm);%repeats the map if there are mulitple cycles per image
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
    varargout{1} = sum(ld,1);
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
    clear datat imap
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
    clear data2 ivec
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
jmax = ngr; 
end

function [data_smoothed] = boxcar_averager(data,reach)
data_smoothed = zeros(size(data));
[~,rows,cols] = size(data(1,:,:));
for i=1+reach:rows-reach
    for j = 1+reach:cols-reach
        for ip = -reach:reach
            for jp = -reach:reach                
                data_smoothed(:,i,j) = data(:,i+ip,j+jp)+data_smoothed(:,i,j);                
            end
        end
    end
end

end

function [new_filename] = remove_sdt(filename)

if length(filename)>3 %fixes filename if sdt is appended to file name
    if strcmp(filename(end-3:end),'.sdt')
        filename=filename(1:end-4);
    end
end
new_filename = filename;

end

function [data128] = make128(data256)
data128 = zeros(4096,128,128);
for i=1:128
    for j=1:128
        i_beg = 2*i-1;
        j_beg = 2*j-1;
        data128(:,i,j) = sum(sum(data256(:,i_beg:i_beg+1,j_beg:j_beg+1),3),2);
    end
end

end
