%Used to Edit/Reset Mlist the MatList

change = 2;
max = 3018;
rmax = 3018;





tempf=load('Z:\bkaye\cluster\matlist.mat','-mat','mlist');
mlist = tempf.mlist;

if change ==1;
clear mlist;
mlist.name = 1:max;
mlist.read(1:rmax) = 1;
mlist.read(rmax+1:max) =0;
save('Z:\bkaye\cluster\matlist.mat','mlist');
end


if change ==2
mlist.name = mlist.name(1:max);
mlist.read = mlist.read(1:max);
save('Z:\bkaye\cluster\matlist.mat','mlist');
end