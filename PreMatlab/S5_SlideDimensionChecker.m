raw_process_dir = 'TonsilProcessing';
group_dir = 'GroupROI';% Name of the folder were organized data will be stored
split_dir = fullfile(data_root, raw_process_dir, group_dir);
mkdir(split_dir); % Makes the directory were the organized data will be stored

mat_dir = 'MatROI'; % Name of the folder with the tissues/patients
slide_list = dir(fullfile(data_root, raw_process_dir, mat_dir)); % list of the slides 
slide_num = length(slide_list) - 2 ; % Length of above list 
for i = 1:slide_num % Loops through each slide
    disp(i + "/" + slide_num);
    p_id = slide_list(i+2).name;
    slide_ROI_dir = fullfile(data_root, raw_process_dir, mat_dir, p_id);
    roi_list = dir(slide_ROI_dir); % list of all the ROIs 
    number_of_roi = length(roi_list) - 2;
    dimensions_list_named = strings(1,3);
    
    % Within each patient folder we will count the total number of files 
    for j = 1:number_of_roi % Loops through the number of ROI in a slide
        roi_id = roi_list(j+2).name; % Gets one ROI at a time to check inside of
        cur_ROI_dir = fullfile(data_root, raw_process_dir, mat_dir, p_id, roi_id);
        image_list = dir(cur_ROI_dir);
        image_list_name = {image_list.name};
        % Creates a list of dimensions with respective ROI 
        load(fullfile(data_root, raw_process_dir, mat_dir, p_id, roi_id,  string(image_list_name(3))), 'stain_img');
        dimensions_list_named(j,1:2) = size(stain_img); 
        dimensions_list_named(j,3) = roi_id;
    end
    
    dimensions_list = dimensions_list_named(:,1:2); % Removes name from dimension list for processing
    sort_categories_list = zeros(1,3); % Empty list to organize categories of dimensions 
    group_list = [0 0]; % Empty lists that hold the types of categories 0 0 is not a category it just simplfies processing
    category_counter = 0; % Counter for category numbers
    
    for a = 1:number_of_roi % Loops through all ROIs and sorts them
        if  ismember(str2double(dimensions_list(a,1:2)),group_list, 'rows') % creates categories for new dimention types
          [match, result]=ismember(group_list,str2double(dimensions_list(a,1:2)),'rows');
          answer = find(result == 1);
          sort_categories_list(a,3) = answer;
          sort_categories_list(a,1:2) = str2double(dimensions_list(a,1:2));
        else % labels ROIs with existing categories so they can be organized
          category_counter = category_counter + 1;
            if category_counter >= 2
              disp("Mismatch Found " + p_id)  
            end 
          group_list(category_counter,:) = str2double(dimensions_list(a,1:2));
          [match, result]=ismember(group_list,str2double(dimensions_list(a,1:2)),'rows');
          answer = find(result == 1);
          sort_categories_list(a,3) = answer;
          sort_categories_list(a,1:2) = str2double(dimensions_list(a,1:2));
        end 
     end 
     disp(p_id + " was checked and mismatches were sorted")
     
    % Uses the number of categories to create new slide folders for each
    % category. 
    p_id_name_list = strings(); 
    for c = 1:category_counter
        p_id_new_name = append(p_id,"-000",int2str(c));
        p_id_name_list(c) = p_id_new_name;
        mkdir (split_dir, p_id_new_name);
    end
    
    % Moves the files from unorganized folder into the new folders created
    % per their dimension category
    for t = 1:number_of_roi
        des_num = sort_categories_list(t,3); % destination number 
        destination_dir = fullfile(data_root,raw_process_dir,group_dir, p_id_name_list(des_num));
        source_dir = fullfile(data_root,raw_process_dir,mat_dir, p_id, roi_list(t+2).name);
        copyfile(source_dir, fullfile(destination_dir, roi_list(t+2).name));
    end 
    
    disp(p_id + " organization completed!")
end
disp("Slide organization completed!")
