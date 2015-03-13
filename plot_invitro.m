    
clear;
imin = 3199;
imax = 3228;
for i = imin:imax
    ind =i-imin+1;
    try
        nstr = strcat('matin',num2str(i),'.mat');
        tempf=load(nstr,'-mat','input');
    catch exception
        nstr = strcat('Z:\bkaye\cluster\matin\matin',num2str(i),'.mat');
        tempf=load(nstr,'-mat','input');
    end
    input = tempf.input(:,:,:);
    dataname = input.dataname;
    %fprintf('%s\n',dataname);
    into(ind) = sum(input(1,1,1).ni./input(1,1,1).sinti);
end
    
int = reshape(into,6,5);
%mint = mean(int,1);





 
intm = mean(int,1);
intv = std(int,0,1);

intmin = max(intm)/2^(length(intm)-1);
intx = intmin*2.^(0:length(intm)-1);
figure(1); clf; hold all;
plot(intx,intm,'bo'); plot(intx,intx,'r'); errorbar(intx,intm,intv);