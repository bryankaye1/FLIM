
function mod_matin(load_nums,varargin)

p = inputParser;
p.addOptional('pass_field',{'NA'}, @iscell); % optional arguments are read in next

p.addRequired('load_nums');
p.addParameter('verbose',1,@(x) x==0 || x==1) %then param pairs.
p.addParameter('set_matnum',0, @isnumeric )

p.addParameter('tbac',-pi, @isnumeric )
p.addParameter('tfw',-pi, @isnumeric )
p.addParameter('w2',-pi, @isnumeric )
p.addParameter('w1',-pi, @isnumeric )
p.addParameter('fracstep',-pi, @isnumeric )
p.addParameter('datahis',-pi, @isnumeric )
p.addParameter('pass_image',-pi)

p.addParameter('w1range',-pi, @isnumeric )


%You cannot enter in param pairs if you did not enter optional arguments.
p.parse(load_nums,varargin{:});
fn_param = p.Results;


%Optional args: verbose, set_matnum
%value pairs: tbac, tfw, w2min,w1max,

%This script is used to change the input parameters of a matin file. Set
%create_new = 1 to create a new matin, or 0 to replace the old matin.

%%load_nums=21465:21501;

%load_beg = 15734;
%load_end = 15789;
set_matnum = fn_param.set_matnum;
verbose = fn_param.verbose;

ind = 0;
input = 1;
for i = fn_param.load_nums%load_beg:load_end
    ind = ind+1;
    [input,~,~] = load_mat_data(i,0);
    if iscell(fn_param.pass_image)
        int_image = fn_param.pass_image;
        seg_results = int_image;
        int_image_flag = 0;
        TEST_FOR_ERROR = int_image{1};
    else
        try
        [int_image,seg_results] = load_int_image(i);
        catch
        end
    end
    
    jmax = input(1,1,1).jmax;
    exptmax = input(1,1,1).exptmax;
    cyclesmax = input(1,1,1).cyclesmax;
    
    for jind = 1:jmax
        for expt = 1:exptmax
            for cindex = 1:cyclesmax
                
                if fn_param.tbac ~=-pi
                    input(cindex,expt,jind).tbac = fn_param.tbac;
                end
                if fn_param.tfw ~=-pi
                    input(cindex,expt,jind).tfw = fn_param.tfw;
                end
                if fn_param.w1 ~=-pi
                    input(cindex,expt,jind).w1min = fn_param.w1;
                    input(cindex,expt,jind).w1max = fn_param.w1;
                end
                
                if fn_param.w2 ~=-pi
                    input(cindex,expt,jind).w2min = fn_param.w2;
                    input(cindex,expt,jind).w2max = fn_param.w2;
                end
                
                if fn_param.w1range ~=-pi
                    input(cindex,expt,jind).w1min = fn_param.w1range(1);
                    input(cindex,expt,jind).w1max = fn_param.w1range(2);
                end
                
                if fn_param.fracstep ~=-pi
                    input(cindex,expt,jind).prstep = fn_param.fracstep;
                    input(cindex,expt,jind).w02step = fn_param.fracstep;
                end
                
                if fn_param.datahis ~=-pi
                    input(cindex,expt,jind).datahis = fn_param.datahis;
                end
                
                if ~strcmp(fn_param.pass_field{1},'NA')
                    input(cindex,expt,jind).([fn_param.pass_field{1}])...
                        = fn_param.pass_field{2};
                end
                
            end
        end
    end
    
    if set_matnum==0
        [MatName,SimName] = write_to_mlist;
        set_matnum=0;
    else
        [MatName,SimName] = write_to_mlist(set_matnum+ind-1);
    end
    if i == load_nums(1)
        newmat_start = MatName(35:strfind(MatName,'.')-1);
    end
    
    save(MatName);
    
    
end
if verbose
    fprintf('matin %s:%s was modified and saved to matin %s:%s\n',...
        num2str(load_nums(1)),num2str(load_nums(end)),...
        newmat_start,MatName(35:strfind(MatName,'.')-1));
end
end













%     if ~int_image_flag
%         if iscell(fn_param.pass_image)
%             save(MatName, 'input', 'int_image','seg_results');
%         else
%             save(MatName, 'input', 'int_image');
%         end
%     else
%         
%     end

