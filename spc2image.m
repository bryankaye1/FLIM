function [int_image] = spc2image(tmini,tmaxi,file_name,pth_sdt,FOV_fn,FOV_pth,tsm,reach)
%Loads spc and returns linearized FLIMage intensity corrected photon counts
%It can return an intensity image if varargin{1} is 'save_int_image'

if tsm ==0
    tsm = 1;
end
for ts=1:tsm
    file_name = remove_sdt(file_name); %removes ".sdt" from end of filename if present
    if tsm ==1
        file_name2 = strcat(file_name,'.sdt');
    else %This subsection padds the suffix '_c0...' with the appropriate amount of zeros.
        %could have used sprintf instead of if statements
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
    
    %datat(:,:,1+128*(ts-1):128*ts) = ld(tmini:tmaxi,:,:);
    %datat = ld(tmini:tmaxi,:,:);
    nph = squeeze(sum(ld(tmini:tmaxi,:,:),1));
    %hislen = tmaxi-tmini+1;
end

%This section loads and makes the intensity-normalization map
if strcmp(FOV_fn,'none')
    imap1 = ones(4096,128,128); %%If you do not input int map, we use a pi map
    imap1 = sum(imap1(tmini:tmaxi,:,:),1); 
else
    FOV_fn = remove_sdt(FOV_fn);
    FOV_fn = [FOV_fn,'.sdt'];
    imap1 = sdt_to_vector(FOV_pth, FOV_fn);
    imap1 = sum(imap1(tmini:tmaxi,:,:),1); 
    imap1 = imgaussfilt(imap1/max(imap1(:)),4);
end 
    imap1 = squeeze(imap1);
    imap = repmat(imap1,1,tsm);%repeats the map if there are mulitple cycles per image
    imap = imap/mean(mean(imap));
    nph = nph./imap;
    
if reach>0
    datat = boxcar_averager(datat,reach);
    nph = boxcar_averager_int(nph,reach); %NEEDS TESTING
end
int_image = nph;
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

function [data_smoothed] = boxcar_averager_int(data,reach)
data_smoothed = zeros(size(data));
[rows,cols] = size(data);
for i=1+reach:rows-reach
    for j = 1+reach:cols-reach
        for ip = -reach:reach
            for jp = -reach:reach                
                data_smoothed(i,j) = data(i+ip,j+jp)+data_smoothed(i,j);                
            end
        end
    end
end
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


% function [new_filename] = remove_sdt(filename)
% 
% if length(filename)>3 %fixes filename if sdt is appended to file name
%     if strcmp(filename(end-3:end),'.sdt')
%         filename=filename(1:end-4);
%     end
% end
% new_filename = filename;
% 
% end


% dataout = pi*ones(hislen,128+(128*tsm-1)*128);
% nphi = pi*ones(128+(128*tsm-1)*128);
% for k = 1:128*tsm
%     for j = 1:128
%         nphi(j+(k-1)*128) = nph(1,j,k)/imap(1,j,k); %total photon number of data, reshaped into line
%         dataout(:,j+(k-1)*128)=datat(:,j,k);
%     end
% end
