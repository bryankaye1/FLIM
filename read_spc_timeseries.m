function [data, timeinterval] = read_spc_timeseries(filename,ndiv)

fileID = fopen(filename);
data1 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
fclose(fileID);
%%
allmicrot = data1{2};
allmacrot = data1{1};
%%
tstep = floor(allmacrot(end)/ndiv);
ti(1) = 0;
for i = 1:ndiv
    
    ti(i+1) = find(allmacrot>=(tstep*i),1);

    microt = allmicrot(ti(i)+1:ti(i+1));
    edges = 0:4096;
    data(i,:)= histc(microt,edges);
    
    timeinterval(i) = (allmacrot(ti(i+1))-allmacrot(ti(i)+1))*25*10^(-9);
    
    
end
end

%%% Code for divided by photon number instead of time %%%
%npg = 100000
%imax = floor(length(allmicrot)/npg);
% for i = 1:imax
%     
%     microt = allmicrot((i-1)*npg+1:i*npg);
%     
%     np = length(microt);
%     edges = 0:4096;
%     data(i,:)= histc(microt,edges);
%     
%     macrostart = 25*10^(-9)*allmacrot((i-1)*npg+1);
%     macroend = 25*10^(-9)*allmacrot(i*npg);
%     
%     timeinterval(i) = macroend - macrostart;
%     timepoints(i) = macroend;
%     
% end








