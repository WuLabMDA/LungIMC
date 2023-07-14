mat_dir = 'MatROI'; % Name of the folder with the tissues/patients

slide_list = dir(fullfile(data_root, raw_process_dir, mat_dir)); % list of the slides 
slide_num = length(slide_list) - 2 ; % Length of above list 
disp('Start dimension check');
for i = 1:slide_num
    disp(i + "/" + slide_num);
    p_id = slide_list(i+2).name;
    slide_ROI_dir = fullfile(data_root, raw_process_dir, mat_dir, p_id);
    roi_list = dir(slide_ROI_dir); % list of all the ROIs 
    number_of_roi = length(roi_list) - 2;
   
    % Within each patient folder we will count the total number of files 
    for j = 1:number_of_roi
        roi_id = roi_list(j+2).name; % Gets one ROI at a time to check inside of
        cur_ROI_dir = fullfile(data_root, raw_process_dir, mat_dir, p_id, roi_id);
        image_list = dir(cur_ROI_dir);
        image_list_name = {image_list.name};
        number_of_images = length(image_list) - 2;
        dimensions_list_named = strings(1,3);
        
       for k = 1:number_of_images
           load(fullfile(data_root, raw_process_dir, mat_dir, p_id, roi_id,  string(image_list_name(k+2))), 'stain_img');
           dimensions_list_named(k,1:2) = size(stain_img); 
           dimensions_list_named(k,3) = string(image_list_name(k+2));
       end
       
       dimensions_list = dimensions_list_named(:,1:2);    
       for a = 1:number_of_images
           if dimensions_list(1,1:2) ~= dimensions_list(a,1:2)
               disp("There is a mismatch of dimensions within " + roi_id + " in " + p_id);
           end 
       end 
    end
end
disp('Done dimension check');