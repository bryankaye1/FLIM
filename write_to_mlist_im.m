function [MatName,SimName] = write_to_mlist_im

clear mlist

tempf=load('Z:\bkaye\cluster\matlist_im.mat','-mat','mlisti');
mlisti = tempf.mlisti;

index = length(mlisti.name)+1;
mlisti.name(index) = index; 
mlisti.read(index) = 0;

save('Z:\bkaye\cluster\matlisti.mat','mlisti');

MatName = strcat('Z:\bkaye\cluster\matin',num2str(index),'.mat');
SimName = strcat('Z:\bkaye\cluster\sim',num2str(index),'.mat');
end