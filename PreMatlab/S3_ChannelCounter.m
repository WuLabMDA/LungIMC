raw_process_dir = 'TonsilProcessing';
mat_dir = 'MatROI'; % Name of the folder with the tissues/patients
channels_expected = 35; % Input the number of channels that we expect per ROI

slide_list = dir(fullfile(data_root, raw_process_dir, mat_dir)); % list of the slides 
slide_num = length(slide_list) - 2; % Length of above list 

incorrect_roi_num = 0; % Counter for the number incorrect
% For loop will gather the total number of patients and enter each patient folder 
for i = 1:slide_num
    disp(i + "/" + slide_num);
    p_id = slide_list(i+2).name;
    slide_ROI_dir = fullfile(data_root, raw_process_dir, mat_dir, p_id);
    roi_list = dir(slide_ROI_dir); % list of all the ROIs 
    number_of_roi = length(roi_list) - 2;
   
    % Within each patient folder we will count the total number of files 
    for j = 1:number_of_roi
        roi_id = roi_list(j+2).name; % Gets one ROI at a time to check inside of
        cur_ROI_dir= fullfile(data_root, raw_process_dir, mat_dir, p_id, roi_id);
        % Change file type in line below as needed
        mat_filelist = dir(fullfile(cur_ROI_dir, '*.mat')); % Gets the list of files 
        if length(mat_filelist) ~= channels_expected % Checks to make sure expected channels are in the folder
            incorrect_roi_num = incorrect_roi_num + 1;
            disp("There is an error at " + roi_id + " in " + p_id);
            break;
            % disp("There are " + tiff_file_length + " files in the ROI")
        end
        % disp(roi_id + " has the correct number of files")
    end
    % disp(p_id + " completed!")
end 

if incorrect_roi_num == 0 % Outputs result statements to summarize
    disp('Everything looks good, all ROIs have the correct number of channels!');
else
    disp("There are " + incorrect_roi_num + " patients with unexpected number of stainings.");
end