%function [hislong,hisshort, data,dshortt] = read_spc_fn_v2(filename,mexp)
clear
%% Get Long-Lifetime Photon Master List

%% Combine Photons to make distributions
type=1;
if type ==1
    clear data
    edges = 0:4095;
    
    for oloop =1:3
        clear tempdlong.dlong tempdshort.dshort data dlong dshort;
        
        if oloop ==1
            tempdlong = load('ll-at');
            dlong = tempdlong.dlong;
            
            tempdshort = load('sl-at');
            dshort = tempdshort.dshort;
        end
        
        if oloop ==2
            tempdlong = load('lm-at');
            dlong = tempdlong.dlong;
            
            tempdshort = load('sm-at');
            dshort = tempdshort.dshort;
        end
        
        if oloop ==3
            tempdlong = load('lh-at');
            dlong = tempdlong.dlong;
            
            tempdshort = load('sh-at');
            dshort = tempdshort.dshort;
        end
        
        yimax = 8;
        repeatmax = 100;
        
        data(yimax,repeatmax).his=pi;
        data(yimax,repeatmax).nplong = pi;
        data(yimax,repeatmax).npshort = pi; 
        data(yimax,repeatmax).hislong = pi;
        data(yimax,repeatmax).hisshort = pi;
        
        tic
        for yi = 1:yimax
            fprintf('yi is %f\n',yi);
            y = 2^(-yi)*10^8;
            for repeat = 1:repeatmax  
                ratio = 2^(-7);
                
                nplong = round(y*(1-ratio)); %number of long photons
                indlong = randi(length(dlong),1,nplong);
                long = dlong(indlong);
                hislong = histc(long, edges);
                
                npshort = round(y*ratio); % Number of short photons
                indshort = randi(length(dshort),1,npshort); %index for random sampling from short dist
                short = dshort(indshort);
                hisshort = histc(short, edges); % checkratio = log2(length(dshort)/length(dlong))
                
                data(yi,repeat).his=hisshort+hislong;
                data(yi,repeat).nplong = nplong;
                data(yi,repeat).npshort = npshort;
                
                data(yi,repeat).hislong = hislong;
                data(yi,repeat).hisshort = hisshort;
                
            end
        end
        toc

        if oloop ==1
            save('cdye_ratio128x_lowInt','data');
        end
        
        if oloop ==2
            save('cdye_ratio128x_medInt','data');
        end
        
        if oloop ==3
            save('cdye_ratio128x_hiInt','data');
        end    
    end

end


if type ==0
    clear data
    edges = 0:4095;
    
    for oloop =2:3
        clear tempdlong.dlong tempdshort.dshort data dlong dshort;
        
        if oloop ==1
            tempdlong = load('ll-at');
            dlong = tempdlong.dlong;
            
            tempdshort = load('sl-at');
            dshort = tempdshort.dshort;
        end
        
        if oloop ==2
            tempdlong = load('lm-at');
            dlong = tempdlong.dlong;
            
            tempdshort = load('sm-at');
            dshort = tempdshort.dshort;
        end
        
        if oloop ==3
            tempdlong = load('lh-at');
            dlong = tempdlong.dlong;
            
            tempdshort = load('sh-at');
            dshort = tempdshort.dshort;
        end
        
        y = 5*10^7;
        %edges = 0:4095; %number of bins
        
        xmax = 18;
        repeatmax = 100;
        
        data(xmax,repeatmax).his=pi;
        data(xmax,repeatmax).nplong = pi;
        data(xmax,repeatmax).npshort = pi;
        
        data(xmax,repeatmax).hislong = pi;
        data(xmax,repeatmax).hisshort = pi;
        
        tic
        for x = 1:xmax
            fprintf('x is %f\n',x);
            for repeat = 1:repeatmax
                
                ratio = 2^(-x);
                
                nplong = round(y*(1-ratio)); %number of long photons
                indlong = randi(length(dlong),1,nplong);
                long = dlong(indlong);
                hislong = histc(long, edges);
                
                npshort = round(y*ratio); % Number of short photons
                indshort = randi(length(dshort),1,npshort); %index for random sampling from short dist
                short = dshort(indshort);
                hisshort = histc(short, edges); % checkratio = log2(length(dshort)/length(dlong))
                
                data(x,repeat).his=hisshort+hislong;
                data(x,repeat).nplong = nplong;
                data(x,repeat).npshort = npshort;
                
                data(x,repeat).hislong = hislong;
                data(x,repeat).hisshort = hisshort;
                
            end
        end
        toc
        
        %matnum = (x-1)*100+repeat;
        % matname = strcat('Z:\bkaye\cluster\data',num2str(matnum),'.mat');
        if oloop ==1
            save('cdata_low_100rep', 'data');
        end
        
        if oloop ==2
            save('cdata_med_100rep', 'data');
        end
        
        if oloop ==3
            save('cdata_hi_100rep', 'data');
        end
        
    end
    toc
end




%filename = 'Z:\bkaye\Bayes_2014\data\2014-10-7\coumarin_1p4e5_550-88nm_950ex_2000sec_m1.asc';
% fileID = fopen('C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-10-7\eryth-1p5e6_550-88nm-950ex-200sec_m1.asc');
% data1 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
% fclose(fileID);
% 
% dshort = data1{2}; %Long lifetime photons (arrival times): dlong = distribution long
% npshort = length(dshort); %number of photons
% toc
% 
% edges = 0:4095;
% hisshort = histc(dshort, edges);
% 
% save('sm-his','hisshort');
% save('sm-at','dshort');
% 
% %% Get Short-Lifetime Photon Master List%%
% tic
% fileID = fopen('C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-10-7\coumarin_1p1e6_550-88nm_950ex_200sec_m1.asc');
% data2 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
% fclose(fileID);
% dlong = data2{2};
% nplong = length(dlong);
% toc
% 
% edges = 0:4095;
% hislong = histc(dlong, edges);
% 
% save('lm-his','hislong');
% 
% save('lm-at','dlong');
