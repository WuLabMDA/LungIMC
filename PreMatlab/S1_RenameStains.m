roi_dir = 'StudySlidesROIs';

patient_list = dir(fullfile(data_root, roi_dir));
num_patients = length(patient_list) - 2;
for pp = 1:num_patients
    disp(pp + "/" + num_patients);
    p_id = patient_list(pp+2).name;
    raw_p_dir = fullfile(data_root, roi_dir, p_id);
    roi_list = dir(raw_p_dir);
    num_rois = length(roi_list) - 2;
    for rr = 1:num_rois
        roi_id = roi_list(rr+2).name;
        stain_list = dir(fullfile(raw_p_dir, roi_id, '*.ome.tiff'));
        stain_num = length(stain_list);
        for ss = 1:stain_num
            stain_name = stain_list(ss).name;
            stain_name = stain_name(1:end-9);
            if strcmp(stain_name, '147Sm_TIGIT')
                src_stain_path = fullfile(raw_p_dir, roi_id, '147Sm_TIGIT.ome.tiff');
                dst_stain_path = fullfile(raw_p_dir, roi_id, '142Nd_TIGIT.ome.tiff');
                movefile(src_stain_path, dst_stain_path);
            end
            if strcmp(stain_name, '176Yb_NaK-ATPase')
                src_stain_path = fullfile(raw_p_dir, roi_id, '176Yb_NaK-ATPase.ome.tiff');
                dst_stain_path = fullfile(raw_p_dir, roi_id, '176Yb_NaKATPase.ome.tiff');
                movefile(src_stain_path, dst_stain_path);
            end

        end
    end
end
