function [set_matnum] = change_lifetimes(MatName_in,w1vec,int_image,varargin)
numvarargs = length(varargin);
optargs = {0};
optargs(1:numvarargs) = varargin;
[specify_matin] = optargs{:};

%This script is used to change the input parameters of a matin file. Set
%create_new = 1 to create a new matin, or 0 to replace the old matin.
ind = 0;
input = 1;
for i = w1vec

    ind = ind+1;
    load(MatName_in);
    
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
    
    if specify_matin==0
    [MatName,SimName] = write_to_mlist; 
    set_matnum=0;
    else
    matstart = MatName_in(35:strfind(MatName_in,'.')-1);
    [MatName,SimName] = write_to_mlist(str2num(matstart)+ind);
    set_matnum = str2num(matstart)+ind+1;
    end
    
    save(MatName, 'input',int_image);
end

fprintf('w1 is %s:%s:%s ending in Matin %s\n',num2str(w1vec(1)),...
    num2str(w1vec(2)-w1vec(1)),num2str(w1vec(end)),MatName(35:strfind(MatName,'.')-1));

end