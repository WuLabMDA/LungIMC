%% Merge 191Ir and 193Ir

group_dir = 'GroupROI';% Name of the folder were organized data will be stored

group_list = dir(fullfile(data_root, raw_process_dir, group_dir)); % list of the data groups
group_num = length(group_list) - 2 ; % Length of above group list 

for ii = 1:group_num % Loops through each slide
    disp("Prepare group " + ii + " in " + group_num);
    g_id = group_list(ii+2).name;
    org_grp_dir = fullfile(data_root, raw_process_dir, group_dir, g_id);
    roi_list = dir(org_grp_dir);
    num_rois = length(roi_list) - 2;   
    for rr = 1:num_rois
        roi_id = roi_list(rr+2).name;
        org_roi_dir = fullfile(org_grp_dir, roi_id);
        % load 191Ir
        Ir191_mat_path = fullfile(org_roi_dir, '191Ir.mat');
        Ir191_struct = load(Ir191_mat_path);
        Ir191_img = Ir191_struct.stain_img;
        
        % load 193Ir
        Ir193_mat_path = fullfile(org_roi_dir, '193Ir.mat');
        Ir193_struct = load(Ir193_mat_path);
        Ir193_img = Ir193_struct.stain_img; 
        % sum 
        stain_img = Ir191_img + Ir193_img;
        save(fullfile(org_roi_dir, 'Ir191_193'), 'stain_img');        

    end
end
disp("Merge 191Ir and 193Ir completed!")