function [set_matnum] = split_input(input_in,int_image,seg_results,set_matnum,dataname,split_matin)

jmax = input_in(1,1,1).jmax;
input_holdon = input_in;

input = 1;
for k = 1:split_matin
    pstart = 1+(k-1)*(jmax);
    pend = k*(jmax);
    input = input_holdon(1,1,round(pstart):round(pend));
    [MatName,matnum] = write_to_mlist(set_matnum);
    if k==1 % This section prints out the matin #s
        matstart = num2str(matnum); %MatName(35:strfind(MatName,'.')-1);
        matend = num2str(str2num(matstart) + split_matin - 1);
        fprintf('%s split_matin: Matin %s:%s\n',dataname,matstart,matend);
        fileID = fopen('matin_prints.txt','at');
        fprintf(fileID,'DN = %s split_matin: Matin %s-%s\n',dataname,matstart,matend);
        fclose(fileID);
        set_matnum = str2num(matstart);
    end
    set_matnum = set_matnum + 1;
    save(MatName, 'input','int_image','seg_results');
end
end