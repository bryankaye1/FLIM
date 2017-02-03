function [input_str,output,varargout] = load_mat_data(input1,varargin)

p = inputParser;
p.addRequired('i',@isscalar);
p.addOptional('verbose',1,@isscalar);
p.addParameter('elapsed',0,@isscalar);
p.addParameter('local',0,@isscalar);
p.addParameter('cluster',0,@isscalar);
p.addParameter('pause_output_DNE',0,@isscalar);
p.addParameter('load_like',0,@isscalar);
p.parse(input1,varargin{:});
fn_inputs = p.Results;

i = num2str(fn_inputs.i);
% numvarargs = length(varargin);
% optargs = {1};
% optargs(1:numvarargs) = varargin;
% [verbose] = optargs{:};

%Varargin(1) = 1 for warning the matout does not exist

try
    if ispc
        if fn_inputs.local
            nstr_output = ['C:\Users\Bryan\Documents\MATLAB\data\matout\matout',i,'.mat'];
            nstr_input = ['Y:\Users\bkaye\cluster\matin\matin',i,'.mat'];
        else
            nstr_output = ['Y:\Users\bkaye\cluster\matout\matout',i,'.mat'];
            nstr_input = ['Y:\Users\bkaye\cluster\matin\matin',i,'.mat'];
        end
        input_dump = load(nstr_input,'-mat','input');
    else
        nstr_output = ['/Users/bryankaye/Documents/MATLAB/data/matout/matout',i,'.mat'];
        
        if fn_inputs.local
            input_dump = load(['/Users/bryankaye/Documents/MATLAB/data/'...
                'matin/matin',i,'.mat'],'-mat','input');
        elseif fn_inputs.cluster
            input_dump = load(['smb://needlemanfs3.rc.fas.harvard.edu/'...
                'needlemanfs3/Users/bkaye/cluster/matin/matin',i,'.mat'],'-mat','input');
        else
            input_dump.input = 'NA';
        end
    end
    
    
    %load output, conditional load likelihood
    if fn_inputs.load_like
        load(nstr_output,'-mat','output','eltime','like_mat'); %%OLD MATINS (I think 24787 to 25057, the files was named like_mat)
        varargout{3} = like_mat;
    else
        load(nstr_output,'-mat','output','eltime');
    end
    
    flag = 0;
    
    if fn_inputs.elapsed
        fprintf('elapsed time for matout%s was %s hours\n',i,num2str(eltime/3600))
    end
    varargout{2} = eltime;
    
catch exception
    if fn_inputs.verbose
        fprintf('matout%s DNE\n',i);
    end
    flag = 1;
    output = pi;
    varargout{2} = 'N/A';
    if ispc
        nstr = ['Y:\Users\bkaye\cluster\matin\matin',i,'.mat'];
        input_dump = load(nstr,'-mat','input');
    else
        if fn_inputs.local
            input_dump = load(['/Users/bryankaye/Documents/MATLAB/data/'...
                'matin/matin',i,'.mat'],'-mat','input');
        end
    end
    
    
    if fn_inputs.pause_output_DNE
        fprintf('matout %s does not exist.',i);
        input('Continue...?');
    end
end

varargout{1} = flag;
input_str = input_dump.input;
end