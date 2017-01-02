clear;
datahis_comb = 0;
ind = 0;

ind_cell{1} = [24787,24789];
ind_cell{2} = [24798,24801,24808,24811];
ind_cell{3} = [24824,24834,24837];
ind_cell{4} = [24719,24721,24730];

for k = 1:length(ind_cell)
    
    for i = ind_cell{k}
        ind = ind +1;
        [input,~,~] = load_mat_data(i,'verbose',0);
        datahis_comb = input(1,1,1).datahis + datahis_comb;
        
        [imagedata,~] = load_int_image(round(i));
        combined_image{ind} = imagedata;
        fprintf('matnum %s added to combined data\n',input.dataname);
    end
    
    mod_matin(i,'verbose',1, 'datahis',datahis_comb, 'pass_image',combined_image);
end
