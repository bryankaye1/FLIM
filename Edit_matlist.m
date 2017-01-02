%Used to Edit/Reset Mlist the MatList

tempf=load('Y:\Users\bkaye\cluster\matlist.mat','-mat','mlist');
%tempf=load('matlist.mat','-mat','mlist');
mlist = tempf.mlist;



 msize = 8937;
 read = msize-1;
% 
% mlist.read = ones(1,msize);
% mlist.read(read+1:end) = 0;
% mlist.name = mlist.name(1:msize);


save('Y:\Users\bkaye\cluster\matlist.mat','mlist');