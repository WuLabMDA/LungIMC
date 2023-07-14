organize_dir = 'MatROI';

stain_selection = {'139La_CD45RO', '141Pr_aSMA', '142Nd_TIGIT', '143Nd_ICOS', '144Nd_HLA-DR', ...
    '145Nd_CD68', '146Nd_MPO', '148Nd_CD11c', '149Sm_CD73', '150Nd_PD-L1', ...
    '151Eu_CD163', '152Sm_GranzymeB', '153Eu_CD11b', '154Sm_CD14', '155Gd_FoxP3', '156Gd_TIM3', ...
    '159Tb_LAG3', '160Gd_CD31', '161Dy_IDO-1', '162Dy_Ki67', '163Dy_VISTA', '164Dy_B2M', '165Ho_PD-1', ...
    '166Er_CD8a', '167Er_CD33', '168Er_B7-H3', '169Tm_CD45', '170Er_CD94', '171Yb_CD19', '172Yb_CD3e', ...
    '173Yb_CD4', '174Yb_CK', '175Lu_CTLA-4', '176Yb_NaKATPase', '191Ir_191Ir', '193Ir_193Ir'};
stain_marker = cell(1, length(stain_selection));
for ss = 1:length(stain_selection)
    stain_str = stain_selection{ss};
    underline_pos = strfind(stain_str, '_');
    stain_marker{ss} = stain_str(underline_pos(1)+1:end);
end

patient_list = dir(fullfile(data_root, raw_roi_dir));
num_patients = length(patient_list) - 2;
for pp = 1:num_patients
    disp(pp + "/" + num_patients);
    p_id = patient_list(pp+2).name; 
    raw_p_dir = fullfile(data_root, raw_roi_dir, p_id);
    roi_list = dir(raw_p_dir);
    num_rois = length(roi_list) - 2;
    for rr = 1:num_rois
        roi_id = roi_list(rr+2).name;
        tif_save_dir = fullfile(data_root, raw_process_dir, organize_dir, p_id, roi_id);
        if ~exist(tif_save_dir, 'dir')
            mkdir(tif_save_dir)
        end  
        stain_list = dir(fullfile(raw_p_dir, roi_id, '*.ome.tiff'));
        stain_num = length(stain_list);
        for ss = 1:stain_num
            stain_name = stain_list(ss).name;
            underline_pos = strfind(stain_name, '_');
            marker_name = stain_name(underline_pos(1)+1:end-9);            
            if any(strcmp(stain_marker, marker_name))
                stain_path = fullfile(raw_p_dir, roi_id, stain_name);
                stain_ome = bfopen(stain_path);
                stain_img = stain_ome{1}{1};
                save(fullfile(tif_save_dir, marker_name), 'stain_img');
            end
        end  
    end
end