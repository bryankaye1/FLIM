function [MatName,SimName] = write_to_mlist

clear mlist

tempf=load('Y:\Users\bkaye\cluster\matlist.mat','-mat','mlist');
mlist = tempf.mlist;

index = length(mlist.name)+1;
mlist.name(index) = index;
mlist.read(index) = 0;

save('Y:\Users\bkaye\cluster\matlist.mat','mlist');

MatName = strcat('Y:\Users\bkaye\cluster\matin\matin',num2str(index),'.mat');
SimName = strcat('Y:\Users\bkaye\cluster\sim\sim',num2str(index),'.mat');
end