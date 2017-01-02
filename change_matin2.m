
%This script is used to change the input parameters of a matin file. Set
%create_new = 1 to create a new matin, or 0 to replace the old matin.

load_nums=21465:21501;

%load_beg = 15734;
%load_end = 15789;
set_matnum = 0;
verbose = 1;

ind = 0;
input = 1;
for i = load_nums%load_beg:load_end
    ind = ind+1;
    [input,~,~] = load_mat_data(i,0);
    [int_image,int_image_flag] = load_int_image(i);
    
    jmax = input(1,1,1).jmax;
    exptmax = input(1,1,1).exptmax;
    cyclesmax = input(1,1,1).cyclesmax;
    
    for jind = 1:jmax
        for expt = 1:exptmax
            for cindex = 1:cyclesmax
                input(cindex,expt,jind).tbac = 0.0;
                input(cindex,expt,jind).tfw = 0.0;
               % input(cindex,expt,jind).w2min = 3.77;
               % input(cindex,expt,jind).w2max = 3.77;
              %  input(cindex,expt,jind).w1min = 1.0;
              %  input(cindex,expt,jind).w1max = 1.0;
                %input(cindex,expt,jind).prstep = .001;
                %input(cindex,expt,jind).w02step = .001;
            end
        end
    end
    
    if set_matnum==0
        [MatName,SimName] = write_to_mlist;
        set_matnum=0;
    else
        [MatName,SimName] = write_to_mlist(set_matnum+ind-1);
    end
    if ~int_image_flag
        save(MatName, 'input', 'int_image');
    else
        save(MatName, 'input');
    end
    if verbose
        fprintf('matin%s is a modified version of matin%s\n',...
            MatName(35:strfind(MatName,'.')-1),num2str(i));
    end
end



