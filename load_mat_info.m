function [input,output] = load_mat_info(i)

try
    if ispc
        nstr = ['Y:\Users\bkaye\cluster\matout\matout',num2str(i),'.mat'];
        load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
        %fprintf('Number of photons is %s\n', num2str((sum(output(1,1,1).datahis))));
        nstr = ['Y:\Users\bkaye\cluster\matin\matin',num2str(i),'.mat'];
        load(nstr,'-mat','input');
    else
        nstr = ['/Volumes/needlemanfs3/Users/bkaye/cluster/matout/matout',num2str(i),'.mat'];
        load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
        %fprintf('Number of photons is %s\n', num2str((sum(output(1,1,1).datahis))));
        nstr = ['/Volumes/needlemanfs3/Users/bkaye/cluster/matin/matin',num2str(i),'.mat'];
        load(nstr,'-mat','input');
    end
catch exception
    fprintf('matout%s DNE',num2str(i));
end
end