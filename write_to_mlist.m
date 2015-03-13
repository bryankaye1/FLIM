function [MatName,SimName] = write_to_mlist

clear mlist

tempf=load('Z:\bkaye\cluster\matlist.mat','-mat','mlist');
mlist = tempf.mlist;

index = length(mlist.name)+1;
mlist.name(index) = index;
mlist.read(index) = 0;

save('Z:\bkaye\cluster\matlist.mat','mlist');

MatName = strcat('Z:\bkaye\cluster\matin\matin',num2str(index),'.mat');
SimName = strcat('Z:\bkaye\cluster\sim\sim',num2str(index),'.mat');
end