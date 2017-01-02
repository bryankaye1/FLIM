function bryans_fn(input1,varargin)
p = inputParser;
p.addRequired('main_num',@isscalar); %required is read is first
p.addOptional('main_text',{'number'}, @iscell); % optional arguments are read in next
p.addParameter('times_said',1,@(x) x>-1) %then param pairs. 
p.addParameter('ts2',0, @(x) x==0 || x==1 )

%You cannot enter in param pairs if you did not enter optional arguments.
p.parse(input1,varargin{:});
inputs = p.Results;

for i = -2:inputs.times_said
fprintf('%1.2f is my favorite %s\n',inputs.main_num,inputs.main_text{1});
end

end