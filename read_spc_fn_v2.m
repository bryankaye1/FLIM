%function [hislong,hisshort, data,dshortt] = read_spc_fn_v2(filename,mexp)

fileID = fopen(filename);
data1 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
fclose(fileID);

dlong = data1{2};
nplong = length(dlong);
edges = 0:4096;
hislong = histc(dlong,edges);

hisshort = 'pass function read_spc_fn a 1';
data = 'pass function read_spc_fn a 1';
dshortt = 'pass function read_spc_fn a 1';

if mexp == 1
    clear data dshortt;
    fileID = fopen('Z:\bkaye\Bayes_2014\data\2014-10-7\coumarin_1p4e5_550-88nm_950ex_2000sec_m1.asc');
    data2 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
    fclose(fileID);
    
    for i =1:1%18
        ratio = 2^(-i);
        
        npshort = round(nplong*ratio);
        dshortt = data2{2};
        dshort = dshortt(1:npshort);
        hisshort = histc(dshort, edges); % checkratio = log2(length(dshort)/length(dlong))
        data(i,:) = hisshort+hislong;
    end
    
    dshortt = data2{2};
    dshort = dshortt(1:npshort);
    hisshort = histc(dshort, edges);
    
end



if mexp == 2
    clear data dshortt;
    fileID = fopen('Z:\bkaye\Bayes_2013\data\2013-10-11\erthy_2e5_p25pc_m1.asc');
    data2 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
    fclose(fileID);
    
    for i =1:1%18
        ratio = 2^(-i);
        
        for j = 1:30
            npshort = round(nplong*ratio);
            dshortt = data2{2};
            indshort = randi(length(dshortt),1,npshort);
            dshort = dshortt(indshort);
            hisshort = histc(dshort, edges); % checkratio = log2(length(dshort)/length(dlong))
            data(i,:,j) = hisshort+hislong;
        end
    end
    
    dshortt = data2{2};
    dshort = dshortt(1:npshort);
    hisshort = histc(dshort, edges);
    
end



%end