%% Max(NaKATPase,B2M)

group_dir = 'GroupROI';% Name of the folder were organized data will be stored
stain_dir = 'SegROI'; % Name of the folder images used for cell segentation.

group_list = dir(fullfile(data_root, raw_process_dir, group_dir)); % list of the data groups
group_num = length(group_list) - 2 ; % Length of above group list 

for ii = 1:group_num % Loops through each slide
    disp("Prepare group " + ii + " in " + group_num);
    g_id = group_list(ii+2).name;
    org_grp_dir = fullfile(data_root, raw_process_dir, group_dir, g_id);
    seg_grp_dir = fullfile(data_root, raw_process_dir, stain_dir, g_id);
    if ~exist(seg_grp_dir, 'dir')
        mkdir(seg_grp_dir)
    end
    roi_list = dir(org_grp_dir);
    num_rois = length(roi_list) - 2;   
    for rr = 1:num_rois
        roi_id = roi_list(rr+2).name;
        org_roi_dir = fullfile(org_grp_dir, roi_id);
        seg_roi_dir = fullfile(seg_grp_dir, roi_id);
        if ~exist(seg_roi_dir, 'dir')
            mkdir(seg_roi_dir)
        end
        % load NaKATPase
        nak_mat_path = fullfile(org_roi_dir, 'NaKATPase.mat');
        nak_struct = load(nak_mat_path);
        nak_img = nak_struct.stain_img;
        % load B2M
        b2m_mat_path = fullfile(org_roi_dir, 'B2M.mat');
        b2m_struct = load(b2m_mat_path);
        b2m_img = b2m_struct.stain_img; 
        % max 
        mem_max = max(nak_img, b2m_img);
        % normalize
        max_val = prctile(mem_max(:), 99);
        min_val = prctile(mem_max(:), 1);
        mem_max(mem_max > max_val) = max_val;
        mem_max(mem_max < min_val) = min_val;
        mem_img = uint8(255.0 * (mem_max - min_val) / (max_val - min_val));
        imwrite(mem_img, fullfile(seg_roi_dir, 'NaK_B2M_Max.tif'));
        
    end
end
disp("Cell segmentation image preparation completed!")
