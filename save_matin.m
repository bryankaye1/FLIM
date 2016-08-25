
function [set_matnum] = save_matin(input_in,int_image,set_matnum,dataname,w1vec)
input = 1;
input = input_in;
[MatName,~] = write_to_mlist(set_matnum);
matnum = MatName(35:strfind(MatName,'.')-1);
fprintf('%s matin %s\n',dataname,matnum);

fileID = fopen('matin_prints.txt','at');
fprintf(fileID,'DN = %s FN = %s\n',dataname,MatName);
fclose(fileID);
save(MatName, 'input','int_image');
if ~isempty(w1vec)
    set_matnum = change_lifetimes(MatName,w1vec,set_matnum,int_image);
elseif set_matnum
    set_matnum = set_matnum + 1;
end
end