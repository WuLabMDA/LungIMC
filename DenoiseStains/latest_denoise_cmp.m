clearvars;

roi_root_dir = 'E:\LungIMCData\LungROIProcessing\Denoise';
raw_dir = fullfile(roi_root_dir, 'RawROIs');
spillover_dir = fullfile(roi_root_dir, 'SpilloverROIs');
denoise_dir = fullfile(roi_root_dir, 'DenoisedROIs-09-22');
cmp_dir = fullfile(roi_root_dir, 'CmpRawSpilloverDenoise-09-22');


stain_str_list = {'aSMA', 'B2M', 'B7_H3', 'CD3e', 'CD4', 'CD8a', 'CD11b', 'CD11c', 'CD14', 'CD19',...
    'CD31', 'CD33', 'CD45', 'CD45RO', 'CD68', 'CD73', 'CD94', 'CD163', 'CK', 'CTLA_4', 'FoxP3', 'GranzymeB',...
    'HLA_DR', 'ICOS', 'IDO_1', 'Ir191', 'Ki67', 'LAG3', 'MPO', 'NaKATPase', 'PD_1', 'PD_L1', 'TIGIT', 'TIM3', 'VISTA'};
stain_pixel_num = [50, 50, 37, 5, 5, 25, 37, 5, 50, 5, ...
    37, 5, 50, 37, 50, 50, 37, 37, 50, 35, 5, 5, ...
    50, 5, 5, 25, 50, 25, 35, 50, 35, 30, 5, 5, 35];
stain_agg_map = containers.Map(stain_str_list, stain_pixel_num);

quantile_val_num = [0.90, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.70, 0.05, ...
    0.80, 0.05, 0.05, 0.05, 0.05, 0.90, 0.85, 0.05, 0.60, 0.05, 0.05, 0.05, ...
    0.80, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.70];
quantile_val_map = containers.Map(stain_str_list, quantile_val_num);

roi_list = dir(spillover_dir);
roi_list = roi_list(3:end);

ro_num = length(roi_list);
for sind = 1:ro_num
% for sind = 420:540
    if mod(sind, 10) == 0
        disp(['Processing ',num2str(sind),'/',num2str(ro_num)])
    end
    roi_name = roi_list(sind).name;
    raw_roi_dir = fullfile(raw_dir, roi_name);
    input_roi_dir = fullfile(spillover_dir, roi_name);
    denoise_roi_dir = fullfile(denoise_dir, roi_name);
    if ~exist(denoise_roi_dir, 'dir')
        mkdir(denoise_roi_dir)
    end
    cmp_roi_dir = fullfile(cmp_dir, roi_name);
    if ~exist(cmp_roi_dir, 'dir')
        mkdir(cmp_roi_dir)
    end    
    
    img_list = dir(input_roi_dir);
    img_list = img_list(3:end);
    for rind = 1:length(img_list)
        stain_fullname = img_list(rind).name;
        stain_name = stain_fullname(1:end-5);
        stain_path = fullfile(input_roi_dir, stain_fullname);
        stain_img = imread(stain_path);
        % remove top strong signals
        stain_max = prctile(stain_img(:), 99.5);
        stain_img(stain_img > stain_max) = stain_max;        
        % smoothing
        denoise_img = medfilt2(stain_img);
        % removal of small regions
        bw_img = imbinarize(denoise_img);
        img_mask = bwareaopen(bw_img, stain_agg_map(stain_name));
        denoise_img(~img_mask) = 0;
        % remove weak signal
        d_thresh = quantile(denoise_img(denoise_img > 0.0), quantile_val_map(stain_name));
        denoise_img(denoise_img < d_thresh) = 0; 
        denoise_path = fullfile(denoise_roi_dir, stain_fullname);
        imwrite(denoise_img, denoise_path);
        
        % prepare compparison image list 
        raw_roi_path = fullfile(raw_roi_dir, stain_fullname);
        raw_img = imread(raw_roi_path);
        spillover_roi_path = fullfile(input_roi_dir, stain_fullname);
        spillover_img = imread(spillover_roi_path);
        high_thresh = quantile(raw_img(raw_img > 0), 0.96);
        % view and save
        montage({raw_img, spillover_img, denoise_img}, 'Size',[1 3], 'DisplayRange', [0, high_thresh]);
        cmp_roi_path = fullfile(cmp_roi_dir, strcat(stain_name, '.png'));
        imwrite(getframe(gcf).cdata, cmp_roi_path);
    end
end