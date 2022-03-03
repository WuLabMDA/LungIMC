raw_dir = 'PROCESSING';
group_dir = 'GroupROI';% Name of the folder were organized data will be stored
stain_dir = 'SegROI'; % Name of the folder images used for cell segentation.

group_list = dir(fullfile(data_root, raw_dir, group_dir)); % list of the data groups
group_num = length(group_list) - 2 ; % Length of above group list 

for ii = 1:group_num % Loops through each slide
    disp("Prepare group " + ii + " in " + group_num);
    g_id = group_list(ii+2).name;
    org_grp_dir = fullfile(data_root, raw_dir, group_dir, g_id);
    seg_grp_dir = fullfile(data_root, raw_dir, stain_dir, g_id);
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
        % save membrane image
        mem_mat_path = fullfile(org_roi_dir, 'NaKATPase.mat');
        load(mem_mat_path, 'stain_img');
        max_val = prctile(stain_img(:), 99);
        min_val = prctile(stain_img(:), 1);
        stain_img(stain_img > max_val) = max_val;
        stain_img(stain_img < min_val) = min_val;
        mem_img = uint8(255.0 * (stain_img - min_val) / (max_val - min_val));
        imwrite(mem_img, fullfile(seg_roi_dir, 'NaKATPase.tif'));
        % save nuclear image
        nuc_mat_path = fullfile(org_roi_dir, '191Ir.mat');
        load(nuc_mat_path, 'stain_img');
        max_val = prctile(stain_img(:), 99);
        min_val = prctile(stain_img(:), 1);
        stain_img(stain_img > max_val) = max_val;
        stain_img(stain_img < min_val) = min_val;
        nuc_img = uint8(255.0 * (stain_img - min_val) / (max_val - min_val));        
        imwrite(nuc_img, fullfile(seg_roi_dir, '191Ir.tif'));
    end
end
disp("Cell segmentation image preparation completed!")