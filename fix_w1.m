function [w1sliced_mats] = fix_w1(first_group,w1_des,last_mat)
% When you have a shortlifetime sweep fixing the dataset(s) and you want to
% only select a certain lifetime. All datasets must have the same w1sweep. 

%first group is a vector of the matins for the first dataset
%w1_des is the desired w1 (must be contained in the w1sweep)
%last_mat the matin of the last w1 of the last dataset
%Enter code here
[~,out_1,flag1] = load_mat_data(first_group(1),1);
[~,out_end,flag2] = load_mat_data(first_group(end),1);
if flag1 || flag2
    error('Some of your MAT files do not exist');
end

w1_start = out_1(1,1,1).w1min;
w1_end = out_end(1,1,1).w1min;
w1vec = w1_start:(w1_end-w1_start)/(length(first_group)-1):w1_end;
w1_des_ind = find(w1vec==w1_des,1);

w1sliced_mats = first_group(1) + w1_des_ind-1: length(first_group):last_mat;




end