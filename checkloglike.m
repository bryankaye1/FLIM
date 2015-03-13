
function checkloglike(loglike)
if sum(sum(sum(sum(sum(isnan(loglike)))))) > 0
    disp('ERROR ERROR ERROR: NaNs in loglike');
end
b = loglike(loglike==-pi);
if b~=0
    disp('ERROR ERROR ERROR: Did not search whole space');
end

end