function [data] = read_imagethr_14ns(tmini,tmaxi,file_name,pth_sdt,thresh,thrt,DesPho)

data=dir([pth_sdt file_name]); %does this line do anything?
block=1;
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);

sild = size(ld); %checks to see if this is an image (3dimensional) or just a FIFO (2 dimensional)
if length(sild) == 3
    if thrt == 0
        ld = squeeze(sum(sum(ld,3),2)); %Averages all the pixels together
        data = ld(tmini:tmaxi);
    else
        ld1 = squeeze(sum(ld(tmini:tmaxi,:,:),1));
        
        figure(13); clf;  ldp1 = ld1; imagesc(ldp1);
        tbefore = ['Before Filter ' file_name]; title(tbefore);
        
        btph = sum(sum(ld1)); btpx = size(ld1,1)*size(ld1,2);
        fprintf('\nb thresh photons %1.0f and pixs %1.0f\n',btph,btpx);
        
        ld2 = ld1;
        if thrt == 1
            atpx = sum(sum(ld2>(max(max(ld2))/thresh)));
            ld2(ld2<(max(max(ld2))/thresh))=0;
            %         else
            %             atpx = sum(sum(ld2>thresh));
            %             ld2(ld2<thresh)=0;
        end
        
        if thrt==2
            thresh = 50;
            atph = sum(sum(ld2));
            atpx = sum(sum(ld2>(max(max(ld2))/thresh)));
            
            while atph > DesPho
                ldtemp = ld2;
                thresh = thresh - 0.01;
                atpx = sum(sum(ldtemp>(max(max(ldtemp))/thresh)));
                ldtemp(ldtemp<(max(max(ldtemp))/thresh))=0;
                atph = sum(sum(ldtemp));
            end
            ld2 = ldtemp;
            
        end
        
        atph = sum(sum(ld2)); fprintf('A thresh photons: %1.0f and pixs: %1.0f\n',atph, atpx); drawnow;
        figure(14); clf; ldp2 = ld2; imagesc(ldp2);
        tafter = ['After Filter ' file_name]; title(tafter);
        
        rld2 = reshape(ld2,128*128,1);
        rld = reshape(ld,4096,128*128);
        %ldt = rld(:,rld2~=0);
        ld3 = sum(rld(:,rld2~=0),2);
        %squeeze(sum(sum(ld(ld2~=0)),3),2); %Averages all the pixels together
        data = ld3(tmini:tmaxi);
    end
else
    data = ld(tmini:tmaxi);
end


data=data';
end