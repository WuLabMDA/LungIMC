clearvars;

% roi_root_dir = 'E:\LungIMCData\LungROIProcessing\Denoise';
% roi_root_dir = 'E:\LungIMCData\HumanWholeIMC\LungROIProcessing\Denoise';
roi_root_dir = 'E:\LungIMCData\HumanSampling35-0\LungROIProcessing\Denoise';
% roi_root_dir = 'E:\LungIMCData\HumanSampling35-1\LungROIProcessing\Denoise';
% roi_root_dir = 'E:\LungIMCData\HumanSampling35-2\LungROIProcessing\Denoise';
% roi_root_dir = 'E:\LungIMCData\HumanSampling35-3\LungROIProcessing\Denoise';
% roi_root_dir = 'E:\LungIMCData\HumanSampling35-4\LungROIProcessing\Denoise';

spillover_dir = fullfile(roi_root_dir, 'SpilloverROIs');
denoise_dir = fullfile(roi_root_dir, 'DenoisedROIs');

stain_str_list = {'aSMA', 'B2M', 'B7_H3', 'CD3e', 'CD4', 'CD8a', 'CD11b', 'CD11c', 'CD14', 'CD19',...
    'CD31', 'CD33', 'CD45', 'CD45RO', 'CD68', 'CD73', 'CD94', 'CD163', 'CK', 'CTLA_4', 'FoxP3', 'GranzymeB',...
    'HLA_DR', 'ICOS', 'IDO_1', 'Ir191', 'Ki67', 'LAG3', 'MPO', 'NaKATPase', 'PD_1', 'PD_L1', 'TIGIT', 'TIM3', 'VISTA'};
stain_pixel_num = [50, 50, 37, 20, 20, 20, 37, 25, 50, 15, ...
    37, 25, 50, 50, 50, 50, 50, 37, 50, 35, 20, 25, ...
    50, 25, 37, 25, 50, 25, 35, 50, 35, 30, 25, 25, 35];
stain_agg_map = containers.Map(stain_str_list, stain_pixel_num);
% quantile_val = 0.05;
quantile_val_num = [0.95, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.85, 0.05, ...
    0.80, 0.05, 0.05, 0.05, 0.05, 0.90, 0.90, 0.05, 0.90, 0.05, 0.05, 0.05, ...
    0.80, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.70];
quantile_val_map = containers.Map(stain_str_list, quantile_val_num);

roi_list = dir(spillover_dir);
roi_list = roi_list(3:end);

for sind = 1:length(roi_list)
    roi_name = roi_list(sind).name;
    input_roi_dir = fullfile(spillover_dir, roi_name);
    denoise_roi_dir = fullfile(denoise_dir, roi_name);
    if ~exist(denoise_roi_dir, 'dir')
        mkdir(denoise_roi_dir)
    end
    img_list = dir(input_roi_dir);
    img_list = img_list(3:end);
    for rind = 1:length(img_list)
        stain_fullname = img_list(rind).name;
        stain_name = stain_fullname(1:end-5);
        stain_path = fullfile(input_roi_dir, stain_fullname);
        stain_img = imread(stain_path);
        % remove top strong signals
        stain_max = prctile(stain_img(:), 99.95);
        stain_img(stain_img > stain_max) = stain_max;        
        % smoothing
        denoise_img = medfilt2(stain_img);
        % removal of small regions
        bw_img = imbinarize(denoise_img);
        img_mask = bwareaopen(bw_img, stain_agg_map(stain_name));
        denoise_img(~img_mask) = 0;
%         % remove weak signal
%         d_thresh = quantile(denoise_img(denoise_img > 0.0), quantile_val_map(stain_name));
%         denoise_img(denoise_img < d_thresh) = 0; 
        denoise_path = fullfile(denoise_roi_dir, stain_fullname);
        imwrite(denoise_img, denoise_path);
    end
end