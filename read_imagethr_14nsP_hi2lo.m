function [data,jmax, intb] = read_imagethr_14nsP_hi2lo(tmini,tmaxi,file_name,pth_sdt,thresh,thrt,DesPho)

data=dir([pth_sdt file_name]); %does this line do anything?
block=1;
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);

sild = size(ld); %checks to see if this is an image (3dimensional) or just a FIFO (2 dimensional)
if length(sild) == 3
    if thrt == 0
        %ld = squeeze(sum(sum(ld,3),2)); %Averages all the pixels together
        datat = ld(tmini:tmaxi,:,:);
        
        for j=1:128
            for k = 1:128
                pmt(j,k) = sum(datat(:,j,k));
            end
        end
        
        the=max(max(pmt))+1;
        g=0;
        
        while the > min(min(pmt))
            g = g+1;
            ths = the-1;
            pgroup = 1;
            while sum(sum(pgroup)) < DesPho && the > min(min(pmt))
                the = the - 1;
                pmigr = pmt<=ths;
                pmile = pmt>=the;
                
                pmi(g,:,:) = pmigr.*pmile;
                pgroup = pmt.*squeeze(pmi(g,:,:));
            end
            thsi(g) = ths;
            thei(g) = the;
        end
        
        for k = 1:g-1
            intb(k) =  (thei(k) + thsi(k))/2;
            
            datag = 0;
            for l=1:128
                for m =1:128
                    if pmi(k,l,m)==0
                    else
                        datag = datat(:,l,m)+datag;
                    end
                end
            end
            dataout(:,k) = datag;
        end
        
        
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


data=dataout';
jmax = g-1;

end