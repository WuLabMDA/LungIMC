clearvars;


roi_root_dir = 'E:\LungIMCData\LungROIProcessing\Denoise';
raw_dir = fullfile(roi_root_dir, 'RawROIs');
spillover_dir = fullfile(roi_root_dir, 'SpilloverROIs');
denoise_dir = fullfile(roi_root_dir, 'DenoisedROIs');
cmp_dir = fullfile(roi_root_dir, 'CmpRawSpilloverDenoise');

roi_list = dir(raw_dir);
roi_list = roi_list(3:end);
for sind = 1:length(roi_list)
    roi_name = roi_list(sind).name;
    raw_roi_dir = fullfile(raw_dir, roi_name);
    spillover_roi_dir = fullfile(spillover_dir, roi_name);
    denoise_roi_dir = fullfile(denoise_dir, roi_name);
    cmp_roi_dir = fullfile(cmp_dir, roi_name);
    if ~exist(cmp_roi_dir, 'dir')
        mkdir(cmp_roi_dir)
    end
    
    img_list = dir(raw_roi_dir);
    img_list = img_list(3:end);
    for rind = 1:length(img_list)
        stain_fullname = img_list(rind).name;
        stain_name = stain_fullname(1:end-5);
        raw_roi_path = fullfile(raw_roi_dir, stain_fullname);
        raw_img = imread(raw_roi_path);
        spillover_roi_path = fullfile(spillover_roi_dir, stain_fullname);
        spillover_img = imread(spillover_roi_path);
        denoised_roi_path = fullfile(denoise_roi_dir, stain_fullname);
        denoised_img = imread(denoised_roi_path);
        cmp_roi_path = fullfile(cmp_roi_dir, strcat(stain_name, '.png'));
        high_thresh = quantile(raw_img(raw_img > 0), 0.96);
        montage({raw_img, spillover_img, denoised_img}, 'Size',[1 3], 'DisplayRange', [0, high_thresh]);
        imwrite(getframe(gcf).cdata, cmp_roi_path);
    end
end