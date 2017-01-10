
function [set_matnum] = save_matin(input_in,int_image,seg_results,set_matnum,dataname,w1vec)
input = 1;
input = input_in;
[MatName,matnum] = write_to_mlist(set_matnum);
matnum = num2str(matnum);
%matnum = MatName(35:strfind(MatName,'.')-1);

fileID = fopen('matin_prints.txt','at');
fprintf(fileID,'DN = %s FN = %s\n',dataname,MatName);
fclose(fileID);
save(MatName, 'input','int_image','seg_results');
if ~isempty(w1vec)
    fprintf('%s (w1sweep %s:%s:%s).  MATIN %s:%2.0f\n',dataname,num2str(w1vec(1)),...
    num2str(w1vec(2)-w1vec(1)),num2str(w1vec(end)), matnum,str2double(matnum)+length(w1vec));
    set_matnum = change_lifetimes(MatName,matnum,w1vec,set_matnum,int_image);
elseif set_matnum
    fprintf('%s matin %s\n',dataname,matnum);
    set_matnum = set_matnum + 1;
else
    fprintf('%s matin %s\n',dataname,matnum);  
    set_matnum = str2num(matnum)+1; %#ok<ST2NM>
end
end


%fprintf('w1 was added from %s:%s:%s ending in Matin %s\n',num2str(w1vec(1)),...
%    num2str(w1vec(2)-w1vec(1)),num2str(w1vec(end)),MatName(35:strfind(MatName,'.')-1));