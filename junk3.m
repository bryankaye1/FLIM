function [] = change_lifetimes(matin_num,w1vec)
%This script is used to change the input parameters of a matin file. Set
%create_new = 1 to create a new matin, or 0 to replace the old matin.

for i = w1vec
    nstr = strcat('Y:\Users\bkaye\cluster\matin\matin',num2str(matin_num),'.mat');
    load(nstr);
    
    jmax = input(1,1,1).jmax;
    exptmax = input(1,1,1).exptmax;
    cyclesmax = input(1,1,1).cyclesmax;
    
    for jind = 1:jmax
        for expt = 1:exptmax
            for cindex = 1:cyclesmax
                input(cindex,expt,jind).w1min = i;%.46;%0.438;
                input(cindex,expt,jind).w1max = i;%3.729;
            end
        end
    end
    
    [MatName,SimName] = write_to_mlist;
    fprintf('w1 is %s Matname = %s\n',num2str(i),MatName);
    save(MatName, 'input');
end
end