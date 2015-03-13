
cyclesmax =  100;
tmini = 1;
tmaxi = 128;

                for cycle = 1:cyclesmax
                    if cycle <10
                        str(cycle).name = strcat('740nm_4sint-10stl-convolaria_nonpzt_c00',num2str(cycle),'.sdt');
                    elseif cycle <100
                        str(cycle).name = strcat('740nm_4sint-10stl-convolaria_nonpzt_c0',num2str(cycle),'.sdt');
                    elseif cycle <999
                        str(cycle).name = strcat('740nm_4sint-10stl-convolaria_nonpzt_c',num2str(cycle),'.sdt');
                    end
                end
           
for cindex = 1:cyclesmax
            if jind ==1
                
                pth_data = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-7-18\';
                dataname = str(cindex).name;
                %dataname = '740nm_4sint-10stl-convolaria_c001.sdt';
                thrt = 0; %0 means do not threshold, 1 means threshold , 2 means treshold for desired photon amount
                imthr = 2;%amount to threshol by. Thresh is max(int)/imthr
                DesPho = 500000;%jind*100000;
                
                [p] = read_imagethr_14ns_debug(tmini,tmaxi,dataname,pth_data,imthr,thrt,DesPho);
            end
            
            
                        intensity(cindex) = sum(p);
end