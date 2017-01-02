
function [error] = checkloglike(loglike)
error = 0;
if sum(sum(sum(sum(sum(isnan(loglike)))))) > 0
    disp('ERROR ERROR ERROR: NaNs in loglike');
    error = 111;
end
b = loglike(loglike==-pi);
if b~=0
    disp('ERROR ERROR ERROR: Did not search whole space');
    error =123;
end

end