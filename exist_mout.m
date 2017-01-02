tempf=load('matlist.mat','-mat','mlist');
mlist = tempf.mlist;
n= 0;
for i=7168:7551
   
    nstr = strcat('Y:\Users\bkaye\cluster\matout\matout',num2str(i),'.mat');
    if exist(nstr)==0
    fprintf('Matout %i DNE\n',i);
    mlist.read(i) =0;
    n = n+1;
    end

end

save('matlist.mat','-mat','mlist'); 
fprintf('n is %f\n', n);