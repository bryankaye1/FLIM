
DesPho = 30;
%datat = ld(tmini:tmaxi,:,:);
        
        for j=1:5
            for k = 1:5
                pmt(j,k) = j*k;
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
            pgroup
            thsi(g) = ths;
            thei(g) = the;
        end
        
%         for k = 1:g-1
%             intb(k) =  (thei(k) + thsi(k))/2;
%             
%             datag = 0;
%             for l=1:5
%                 for m =1:5
%                     if pmi(k,l,m)==0
%                     else
%                         datag = datat(:,l,m)+datag;
%                     end
%                 end
%             end
%             dataout(:,k) = datag;
%         end