%% Channel detection script - Ruaridh Clark
% This script detects tidal channel paths from a SAR image with reference 
% to previous imagery and detected channel path.

clearvars

year = '2021';                              % Define year folder

folder_img = strcat('Images/',year,'/');    % folder containing images

[files,check,old_path,chnl_means,scale] = setup(folder_img); % Setup function

folder_sv = strcat('Saved/',year,'/');      % Set the output folder

save_output = 0;                            % Set a flag to save the output (0 = no, 1 = yes)

first = 1;                                  % Set the index of the first image to be processed
for i = first : length(files)               % assuming chronological order of files
    
    [save_img,nm,lat,lon] = process_image(i,folder_img,files);  % open image file

    if i == first   
        % Define first path
        user_input = 0;
        [PDImg,new_path,cat_all,takePath,ref] = first_path_def(user_input,save_img,folder_sv,nm,scale);
        first_path=[];
    else            
        % Detect the channel in all conditions using memory of previous path
        [PDImg,new_path,first_path,cat_all,takePath] = route_all_conditions(save_img,prev_img,old_path,cat_all,scale,ref);
    end
    
    [spaced_path] = downsampling_points(new_path,scale);                % Space out the waypoints on the detected channel path
    
    [save_img,file] = create_figure(nm,save_img,folder_img);            % Create figure

    if save_output==1 && (i > first || user_input == 1)
        display_path(i,save_img,new_path,nm)                            % Display and save channel path
        export_coordinates(spaced_path,lat,lon,nm)                      % Save coordinates to file
        mean_channel = mean_channel_intensity(spaced_path,save_img);    % Calculate the mean SAR intensity over the detected channel path
        saving_path(nm,file,new_path,mean_channel,cat_all,first_path,folder_sv)     % Save the detected path to a file
    end
    
    old_path = fix_spaced_path(i,new_path,old_path,ref(:,2),ref(:,1),takePath,first);   % If needed, extend spaced_path with previous path points

    prev_img = PDImg;                                                   % Set PDImg as previous image
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [files,check,old_path,chnl_means,scale] = setup(foldernm)
    %   Setup function
    %
    %     Inputs
    %     - foldernm: The name of the folder containing the images.
    
    %     Outputs
    %     - files: The list of all files in the specified folder with ".tif" extension.
    %     - check: A variable used to check if certain conditions are met.
    %     - old_path: The variable used to store the old path.
    %     - chnl_means: A matrix used to store the mean channel intensities for each image.
    %     - scale: The scaling factor used for image processing.

    files = dir(strcat(foldernm,'*.tif'));
    
    check = 0;
    old_path = [];
    chnl_means = zeros(length(files),1);
    scale = 30;
end

function [save_img,nm,lat,lon] = process_image(i,foldernm,files)
    %   Function to clear previous geotiff image data before opening next  image
    %
    %     Inputs
    %     - i: The index of the current image.
    %     - foldernm: The name of the folder containing the images.
    %     - files: The list of all files in the specified folder with ".tif" extension.
    
    %     Outputs
    %     - save_img: The processed image.
    %     - nm: The name of the current image.
    %     - lat: The latitude of the image pixels.
    %     - lon: The longitude of the image pixels.

    nm = files(i).name; 
    nm = erase(nm,'.tif');

    file = strcat(foldernm,nm,".tif");

    clear geotiffinfo
    [save_img,lat,lon] = input_process_image(file);
    display(nm)
end

function [img,lat,lon] = input_process_image(file)
    %   Function to process geotiff input
    %
    %     Inputs
    %     - file: The name of the image file.
    %
    %     Outputs
    %     - img: The processed image.
    %     - lat: The latitude of the image pixels.
    %     - lon: The longitude of the image pixels.

    [dB,R] = readgeoraster(file,'OutputType','double');
    info = geotiffinfo(file);   % Read the metadata of the GeoTIFF file
    height = info.Height;       % Integer indicating the height of the image in pixels
    width = info.Width;         % Integer indicating the width of the image in pixels
    [cols,rows] = meshgrid(1:width,1:height);   % Create two 2D arrays (cols and rows)
    [x,y] = intrinsicToWorld(R,cols,rows);    % Convert the pixel coordinates to (x,y) using the spatial referencing object R.
    [lat,lon] = projinv(info, x,y);    % Convert (x,y) to (latitude, longitude) using the projection information

    % Clear variables
    vars = {'R','info','height','width','cols','rows','x','y'};
    clear(vars{:});
    clear vars;

    % Normalise the dB scale data
    timg = rescale(dB,0,255,'InputMin',-40,'InputMax',0);
    timg = timg(1:end,1:end);
    img=timg;

    % Create matrix where each element represents the pixel row number
    rowsMat = zeros(size(img,1),size(img,2));
    for i = 1 : size(img,1)
        rowsMat(i,:)=i.*ones(1,size(img,2));
    end
    
    tmp = img;
    img(isnan(img))=255;    % Replace all NaN values
    img=medfilt2(img);      % Apply a median filter
    img(isnan(tmp))=NaN;    % Restore the NaN values
end

function [PDImg,new_path,cat_all,takePath,ref] = first_path_def(user_input,save_img,folder_sv,nm,scale)
    % This function defines the first path to be used in the path finding algorithm,
    % either by asking the user to select a starting point (user_input = 1) or by
    % loading a previously saved path (user_input = 0).
    %
    % Inputs
    % - user_input: 1 to define a new path, 0 to load a previously saved path.
    % - save_img: image.
    % - folder_sv: folder where the saved path is stored.
    % - nm: name of the file.
    % - scale: the scale factor for downsampling.
    %
    % Returns:
    % - img: output image.
    % - new_path: set of coordinates defining new path.
    % - cat_all: historic sand/mud masks.
    % - takePath: number of points missing intensity values.
    % - ref: reference coordinates.

    if user_input == 1          % If user_input is 1, ask the user to select a starting point.
        [x,y] = user_defines_channel(save_img);
        ref = [y,x];
        [old_path,~] = find_reference_path(save_img,x,y,0);
        [PDImg,new_path,~,cat_all,takePath] = route_all_conditions(save_img,ones(3422,1555)*255,old_path,[],scale,ref);
    elseif user_input == 0      % If user_input is 0, load a previously saved path.
        load(strcat(folder_sv,nm,'.mat'))
        takePath=1;
        
        % Set reference start and end of path
        if new_path(1,1)>new_path(end,1)
            new_path=flip(new_path);
        end
        x_def = [new_path(1,2),new_path(end,2)];
        y_def = [new_path(1,1),new_path(end,1)];
        ref = [y_def',x_def'];


        % Create path difference image
        tmp_img = save_img;
        tmp_img(isnan(save_img))=255;        % Replace NaN values with 255
        [smooth_img] = smooth_image(tmp_img,9); % smooth image
        Ikeep2=[1:round(length(new_path)/100):length(new_path),length(new_path)];
        path = new_path(Ikeep2,:);
        idx = sub2ind(size(smooth_img), path(:,1),path(:,2));
        [PDImg,~] = route_diff_img(path,smooth_img(idx),save_img);
    end
end

function [PDImg,new_path,first_path,cat_all,takePath] = route_all_conditions(save_img,prev_img,old_path,cat_all,scale,ref)
    % Function to detect channel paths in all conditions
    %
    %   Inputs:
    %   - save_img: initial image.
    %   - prev_img: previous image.
    %   - old_path: previous path.
    %   - cat_all: historic sand/mud masks.
    %   - scale: scaling factor for downsampling.
    %   - ref: reference coordinates.
    %
    %   Outputs:
    %   - PDImg: Path difference image
    %   - new_path: new path.
    %   - first_path: first path found.
    %   - cat_all: historic sand/mud masks.
    %   - takePath: number of points missing intensity values.

    if ~isempty(cat_all)
        save_img = adjust_matrix_size(cat_all(:,:,1),save_img); % Adjust save_img size to match previous images
    end

    [first_path,~,takePath,PDImg] = detect_channel_path(save_img,prev_img,old_path,1,ref,scale); % Find channel without mask

    ref = update_ref(ref,takePath,first_path);                             % Update ref. path end point

    [img_mask,cat_all] = cat_filter(PDImg,prev_img,save_img,cat_all);      % Add historic sand/mud pixel masks

    [tmp_path,~,~,~] = detect_channel_path(img_mask,prev_img,first_path,2,ref,scale); % Find channel with mask
    
    [new_path,~] = correct_channel_deviations(tmp_path,img_mask,PDImg);    % Correct channel where it traverses high intensity pixels

    plot_channels(PDImg,first_path,tmp_path,new_path)           % Plot channels

    prev_img = adjust_matrix_size(PDImg,prev_img);              % Replace missing pixels with previous image pixels
    PDImg(PDImg==255)=prev_img(PDImg==255);
end
        
function [new_path,img,takePath,img_update] = detect_channel_path(img,prev_img,old_path,attmpt,ref,scale)
    % This function detects a channel path in an image using a reference channel 
    % and previous images.
    %
    % Inputs
    % - img: current image.
    % - prev_img: previous image.
    % - old_path: old path coordinates to be updated.
    % - attmpt: an integer indicating the number of attempts to detect the channel.
    % - ref: reference start/end path coordinates.
    % - scale: the scale factor to be used for downsampling the old path.
    %
    % Outputs
    % - new_path: created path.
    % - img: updated image.
    % - takePath: number of points missing intensity values.
    % - img_update: updated image after masking and downsampling.

    idx = sub2ind(size(img),old_path(:,1),old_path(:,2));
    takePath = sum(isnan(img(idx)));                        % Count points missing intensity values

    [old_path] = downsampling_points(old_path,scale);       % Downsample the old path to reduce the number of points.
    idx = sub2ind(size(img),old_path(:,1),old_path(:,2));

    if attmpt == 1      % On the first attempt, compute the difference between the image and the path
        [img,old_path] = img_diff_from_path(img,idx,old_path,takePath);     % Path diff. image
        [vs_img] = define_vs_img(img,prev_img);                             % Compare with prev_img
    elseif attmpt == 2  % On the second attempt, remove all points that have maximum intensity.
        old_path=old_path(img(idx)<255,:);
        vs_img = zeros(size(img,1),size(img,2));
    end

    [img_update,spaced_init] = mask_img_gradient_descend(img,old_path,vs_img); % Apply a mask based on vs_image and downsample path.

    [new_path] = build_newpath(spaced_init,1,img_update,ref);   % Build the new path

    [new_path] = check_for_flip(new_path,ref(:,1),ref(:,2));    % Check for flip in the new path and correct it if necessary.
end

function [new_path,img_mask] = correct_channel_deviations(new_path,img_mask,img)
    % Function corrects channel deviations by removing segments overlapping 
    % high intensity pixels and adding new segments to the path in their place
    %
    % Input:
    %   - new_path: the original path of the channel
    %   - img_mask: grayscale image with the historic sand/mud mask applied
    %   - img: original grayscale image
    %
    % Output:
    %   - new_path: corrected path following the tidal channel
    %   - img_mask: updated image image with the historic sand/mud mask applied

    img_mask(img_mask<255)=img(img_mask<255);

    idx = sub2ind(size(img_mask),new_path(:,1),new_path(:,2));
    [movmn] = movmean(img_mask(idx),round(size(new_path,1)/20));    % moving average (quarter path length)
    ids=find(movmn>25);             % Find the indices where the moving average exceeds a threshold

    if ~isempty(ids)                % check for points exceeding moving average
        d_ids = diff(ids);          % difference between adjacent values
        
        % find segment start/end points
        ids_gt1 = find(d_ids>1);
        ids_end = [ids(ids_gt1);ids(end)];
        ids_start = [ids(1);ids(ids_gt1+1)];
        
        keep=[];
        for j = 1 : length(ids_start) 
            if ~isequal(ids_start(j),ids_end(j))
                keep = [keep,j];
            end
        end
        ids_start = ids_start(keep,:);
        ids_end = ids_end(keep,:);

        for j = 1 : length(ids_start) 
            % Build new path for segment to be removed
            [add_path] = build_newpath(new_path([ids_start(j),ids_end(j)],:),1,img,new_path([ids_start(j),ids_end(j)],:));
            
            % Replace segment with new path
            if j == 1 && j == length(ids_start)
                tmp_path = [new_path(1:ids_start(j)-1,:);add_path;new_path(ids_end(j)+1:end,:)];
            elseif j == 1
                tmp_path = [new_path(1:ids_start(j)-1,:);add_path];
            elseif j == length(ids_start)
                tmp_path = [tmp_path;new_path(ids_end(j-1)+1:ids_start(j),:);add_path;new_path(ids_end(j)+1:end,:)];
            else
                tmp_path = [tmp_path;new_path(ids_end(j-1)+1:ids_start(j),:);add_path];
            end
        end

        new_path = remove_loops(tmp_path);                      % remove loops in path
    end
    
end

function [] = plot_channels(img,first_path,tmp_path,new_path)
    % Function to display the first, second, and final path on path diff. image
    figure;imshow(rescale(img,0,1))  
    hold on                                                     
    scatter(first_path(:,2),first_path(:,1),'m.')               
    hold on                                                     
    scatter(tmp_path(:,2),tmp_path(:,1),'r.')                   
    hold on                                                                 
    scatter(new_path(:,2),new_path(:,1),'y.')  
end

function [img,cat_all] = cat_filter(img,prev_img,save_img,cat_all)
    % Function to create and store historic sand/mud mask
    %
    % Inputs:
    %   - img: Input image.
    %   - prev_img: Previous image.
    %   - save_img: Image to be saved.
    %   - cat_all: Previous sand/mud masks.
    % Outputs:
    %   - img: Filtered image.
    %   - cat_all: Updated sand/mud masks.

    img_tmp = img;
    prev_img = adjust_matrix_size(img,prev_img);
    img_tmp(isnan(save_img))=255;
    img_tmp(img_tmp==255)=prev_img(img_tmp==255);       % Replace NaN values with previous image pixels

    [cat] = create_cat(img_tmp, 40);                    % Create sand/mud mask with threshold at 50
    if ~isempty(cat_all)
        cat = adjust_matrix_size(cat_all(:,:,1),cat);   % Adjust mask to previous
        cat_all(:,:,size(cat_all,3)+1) = cat;           % Add to mask collection
    else
        cat_all(:,:,1) = cat;                           % Start mask collection
    end
    
    
    if size(cat_all,3)<5    % If the number of previous masks is less than 5, average all the masks
        filter = (sum(cat_all,3))./(size(cat_all,3));   
    else
        filter = (sum(cat_all(:,:,end-4:end),3))./(size(cat_all(:,:,end-4:end),3));  
    end
    img(filter>=0.5)=255;   % 50% occurence of sand/mud is the threshold
end

function [cat] = create_cat(img, thresh)
    % Function to create historic sand/mud mask where thresh sets the
    % intensity threshold for difference from tidal channel pixel intensity
    % for defining sand/mud

    cat = img;
    cat(img<=thresh)=0;
    cat(img>thresh)=1;
end

function [PDImg,path] = img_diff_from_path(img,idx,old_path,takePath)
    %% Create image relative to path
    % Function creates path difference image relative to channel path
    %
    % Inputs:
    %   img: input image.
    %   idx: indices of old_path in img.
    %   old_path: path to be used as reference.
    %   takePath: number of points missing intensity values.
    %
    % Outputs:
    %   img: updated image.
    %   path: downsampled path with points removed over high intensity regions.

    old_path(isnan(img(idx)),:)=[]; % Remove points from old_path that are NaN

    tmp_img = img;
    tmp_img(isnan(img))=255;        % Replace NaN values with 255
    [smooth_img] = smooth_image(tmp_img,9); % smooth image

    idx = sub2ind(size(tmp_img), old_path(:,1),old_path(:,2));
    riv_rollmn = [old_path(:,1),old_path(:,2),tmp_img(idx)];    % path and associated intensity
    
    [check,~] = check_SAR(smooth_img,riv_rollmn);  % Compare peak low intensity with river channel median
    
    if check == 1       % if channel is NOT lowest intensity feature
        path = river_path_finder(old_path,smooth_img,median(smooth_img(idx)));
    elseif check == 0   % if channel is lowest intensity feature
        path = river_path_finder(old_path,smooth_img,0);
    end

    [path,idx] = point_joiner(path,smooth_img);     % Linearly join paths to connect reference path

    path((tmp_img(idx)==255),:)=[]; % Remove points that correspond to NaN values in tmp_img
    idx(tmp_img(idx)==255)=[];
    
    % Compute custom moving average of smoothed image
    [movmn] = movmean_pre_stdrmv(smooth_img(idx),round(size(path,1)/4));    % rolling average, only including values within 1 std dev

    diff=smooth_img(idx)-movmn;     % Difference between smoothed image and moving average
    diff(1)=0;                      % keep first and last points
    if takePath>50
        diff(end)=0;
    end
    stddev = std(smooth_img(idx));

    % Find indices within 1 std deviation
    Ikeep = find(abs(diff)<stddev);
    Ikeep = unique([1;Ikeep;length(diff)],'stable');
    path = path(Ikeep,:);
    Lidx = idx(Ikeep);
    movmn = movmn(Ikeep);

    % Downsample
    Ikeep2=1:round(length(movmn)/100):length(movmn);
    Ikeep2=[Ikeep2,length(movmn)];
    path = path(Ikeep2,:);

    [PDImg,~] = route_diff_img(path,smooth_img(Lidx(Ikeep2)),img);
end

function [path,idx] = point_joiner(old_path,smooth_img)
    % This function takes a reference path old_path and the smoothed image 
    % smooth_img, and joins the points along the path to create a 
    % continuous linear path, returning the coordinates of the joined path 
    % and their corresponding linear indices in the smoothed image.
    % 
    % Inputs:
    % - old_path: a matrix of size Nx2 containing the row and column 
    % coordinates of a reference path.
    % - smooth_img: smoothed image.
    %
    % Outputs:
    % - path: coordinates of the joined path.
    % - idx: linear indices of the joined path.
    
    x_all = [];
    y_all=[];        
    %% Doesn't include path to ref end point
    for ii = 1 : size(old_path,1)-2
        % calculate the gradient and y-intercept
        grad=(old_path(ii+1,1)-old_path(ii,1))/(old_path(ii+1,2)-old_path(ii,2));
        c = old_path(ii,1)-grad*(old_path(ii,2));
        if abs(grad)<1  % if the gradient is less than 1, create a sequence of x-coordinates
            if old_path(ii,2)<old_path(ii+1,2)
                shft = 1;
            else
                shft = -1;
            end
            x = old_path(ii,2):shft:old_path(ii+1,2);
            y = round(grad.*x+c);
        else            % create a linear path of x and y points
            if old_path(ii,1)<old_path(ii+1,1)
                shft = 1;
            else
                shft = -1;
            end
            y = old_path(ii,1):shft:old_path(ii+1,1);
            x = round((y-c)/grad);
        end

        y_all = [y_all;y'];
        x_all = [x_all;x'];
    end
    % add the last point of the reference path to the joined path
    y_all = [y_all;old_path(end,1)];
    x_all = [x_all;old_path(end,2)];

    % remove any NaN values in the joined path
    Irmv = [find(isnan(x_all)),find(isnan(y_all))];
    x_all(Irmv)=[];
    y_all(Irmv)=[];

    path = [y_all,x_all];
    idx = sub2ind(size(smooth_img),path(:,1),path(:,2));

    % apply gradient descent 
    path = river_path_finder(path,smooth_img,median(smooth_img(idx)));
    idx = sub2ind(size(smooth_img),path(:,1),path(:,2));
end

function [PDImg,riv_rollmn] = route_diff_img(path,pxls_in,img)
    % Function to associate path points with moving mean intensities and
    % then create path difference image

    [movmn] = movmean(pxls_in,round(size(path,1)/4));   % Moving mean
    riv_rollmn = [path(:,1), path(:,2), movmn];         

    [PDImg] = river_ref_diff(img,riv_rollmn);    % Find difference image
    PDImg= rescale(PDImg,1,255); 
end

function [movmn] = movmean_pre_stdrmv(vector,window)
    % Function calculates the moving mean of a vector while removing points 
    % outside of the specified standard deviation threshold.

    for i = 1 : length(vector)
        if i>window
            vec = vector(1:i-1);
            stddev = std(vec(~isnan(vec))-movmn(~isnan(vec)));
            if abs(vector(i)-movmn(i-1))>stddev     % Remove points outside of threshold
                vector(i)=NaN;
            end
            vec = vector(i-window:i);
            movmn(i,1) = mean(vec(~isnan(vec)));    % Find mean using only non-NaN values
        else
            movmn(i,1) = mean(vector(1:i));
        end
    end
end

function [dist] = dist_compare(vecA,vecB)
    % Function to calculate Euclidean distance

    dist = sqrt((vecA(1)-vecB(1))^2+(vecA(2)-vecB(2))^2);
end

function [ref] = update_ref(ref,takePath,first_path)
    % Function to create new reference point for end of path, 
    % if the full tidal channel is not visible in the image 

    if takePath > 50
        ref(2,:) = first_path(end,:);   % set new reference for the route end
    end
end

function [x,y] = user_defines_channel(img)
    % User defines start and end of channel path

    figure;
    pause(1)
    imshow(rescale(img,0,1));

    [x,y]=ginput(2);    % user input
    x = round(x);
    y = round(y);
end

function [path,range] = find_reference_path(img,x,y,river_ref)
    % Given an image `img`, start and end points `x` and `y`, and a reference
    % river image `river_ref`, this function finds the optimal path between
    % the start and end points that passes through the river, by gradually
    % increasing the river threshold `range` until a successful path is found.
    %
    % Inputs:
    %     img: input image.
    %     x: x-coordinates of the start and end points.
    %     y: y-coordinates of the start and end points.
    %     river_ref: reference image intensity to be used for thresholding.
    %
    % Outputs:
    %     path: reference path between the start and end points
    %     range: threshold used for reference path

    range=1; 
    success=0;
    start = [round(y(1)),round(x(1))];
    finish = [round(y(2)),round(x(2))];
    while success == 0
        % Threshold to create binary matrix of traversable pixels
        [Temp] = image_threshrange(img, river_ref, range);
        M=[];
        % Solve the maze using the thresholded image, start and end points
        [path,success,~] = solve_maze_iter(Temp,start,finish,M,1);
        if success == 1
            break
        end
        range = range + 5;
    end
end

function current_img = adjust_matrix_size(ref_img,current_img)
    % Function adjusts the size of current_img to match the size of ref_img, 
    % by either cropping or padding with zeros. 
    %
    % Inputs:
    %     current_img: current image.
    %     ref_img: reference image.
    %
    % Outputs:
    %     current_img: update current image.

    if size(ref_img,1)<size(current_img,1)
        current_img = current_img(1:size(ref_img,1),:);
    elseif size(ref_img,1)>size(current_img,1)
        add = ref_img(end-(size(ref_img,1)-size(current_img,1))+1:end,:);
        if size(add,2)<size(current_img,2)
            add = [add,zeros(size(add,1),size(current_img,2)-size(add,2))];
        elseif size(add,2)>size(current_img,2)
            add = add(1,1:size(current_img,2));
        end
        current_img = [current_img;add];
    end

    if size(ref_img,2)<size(current_img,2)
        current_img = current_img(:,1:size(ref_img,2));
    elseif size(ref_img,2)>size(current_img,2)
        current_img = [current_img,ref_img(1:size(current_img,1),size(ref_img,2)-size(current_img,2))];
    end
end

function [optns] = reachable_points(bin_mat,init_path)
    % Function to remove points that lie on 0 entries in the binary matrix 
    % of traversable points

    idx = sub2ind(size(bin_mat), init_path(:,1),init_path(:,2));
    optns = init_path(bin_mat(idx)>0,:);    % reachable waypoint
end

function [init_path] = rmv_frm_path(init_path,mid_start,starters,mid_finish,enders)
    % Function removes points from an initial path that are between the midpoint-start
    % and -end points. Also removes the midpoint start and end points if they
    % have been used as both a start and end point for a new path segment. 
    %
    % Inputs:
    %   - init_path: initial path.
    %   - mid_start: [x,y], midpoint start point.
    %   - starters: used starting points for path segments.
    %   - mid_finish: [x,y], midpoint end point.
    %   - enders: used ending points for path segments.
    %
    % Outputs:
    %   - init_path: updated path.

    % Check if point is both start and end of path segment
    Istart = find(init_path(:,1)==mid_start(1) & init_path(:,2)==mid_start(2));
    Iend = find(init_path(:,1)==mid_finish(1) & init_path(:,2)==mid_finish(2));

    % Remove all points between the start and end of the path. Also remove
    % the start and end points if they have been both a start and end point.
    if ~isempty(Istart) && ~isempty(Iend)
        if (max(sum(ismember(starters,mid_start),2))==2 && max(sum(ismember(enders,mid_start),2)==2)) || (Istart == 1 && size(init_path,1)>2)
            start_rmv = Istart;         % 
        else
            start_rmv = Istart+1;
        end
        if (max(sum(ismember(starters,mid_finish),2))==2 && max(sum(ismember(enders,mid_finish),2))==2) || (Iend == length(init_path) && size(init_path,1)>2)
            end_rmv = Iend;
        else
            end_rmv = Iend-1;
        end
    
        if start_rmv<=end_rmv
            init_path(start_rmv:end_rmv,:)=[];
        end
    end
end

function [vs_img] = define_vs_img(img,prev_img)
    % Function to create an image that is the difference between the img
    % and prev_img

    prev_img = adjust_matrix_size(img,prev_img);
    img(isnan(img))=255;
    prev_img(isnan(prev_img))=255;
    vs_img = img-prev_img;
end

function [img,spaced_init] = mask_img_gradient_descend(img,old_path,vs_img)
    % Function performs gradient descent and downsampling on old_path, as
    % well as applying a mask to the pixel with the greatest intensity
    % increase from the previous image

    tmp_img = img;
    img(isnan(img))=255;
    
    old_path = check_path_in_image(img,old_path);

    [~,refPixels] = find_river_ref(tmp_img,old_path);

    spaced_init = old_path(~isnan(refPixels),:);
    
    [smooth_img] = smooth_image(img,9);    % smooth image
    spaced_init = unique(river_path_finder(spaced_init,smooth_img,0),'stable','rows');  % % Gradient descent into river pixels
    
    img(vs_img>75) = 255;                   % Update image with recent bright pixel mask
end

function [old_path] = check_path_in_image(img,old_path)
    % Function ensures old_path within the boundaries of the image img

    A=old_path(:,1);
    B=old_path(:,2);
    A(A>size(img,1))=size(img,1);
    B(B>size(img,2))=size(img,2);
    old_path = [A,B];
end

function [river_ref,refPixels] = find_river_ref(img,old_path)
    % Function identifies image intensities for old_path points and 
    % calculates river_ref as the median pixel intensity

    I = sub2ind(size(img),old_path(:,1),old_path(:,2));
    refPixels = img(I);
    river_ref = median(refPixels(~isnan(refPixels)));
end

function [new_path] = remove_loops(new_path)
    % Removes loops in the path represented by new_path.
    % Input: 
    %   - new_path - input path
    % Output:
    %   - new_path - the modified path after removing loops

    % Create a binary image of the path
    path_img = zeros(max(new_path(:,1))+1,max(new_path(:,2))+1);
    I = sub2ind(size(path_img),new_path(:,1),new_path(:,2));
    path_img(I)=1;
    
    % Find a direct route through the path, removing loops
    [temp_path,success,~,~] = solve_maze_iter(path_img,new_path(1,:),new_path(end,:),[],2); % 1 for initial, 2 for refined?
    
    % If the algorithm was successful, update new_path
    if success==1
        new_path = temp_path;
    end
end

function [seg_new,init_path,enders,starters] = connect_path_par(Temp,optns,seg_new,init_path,enders,starters)
    %% Parallel loop through all available starting points to find route to elligible end point
    % Check for possible route from closest route points to furthest, until successful

    % Function tries to connect two unconnected segments of the path. 
    % It loops through all available starting points to find the route to 
    % an eligible end point, and then saves the details of the new segment. 
    %
    % Inputs:
    % - Temp: A binary image representing the mask of the river.
    % - optns: A matrix containing the coordinates of candidate start and end points for path segments.
    % - seg_new: A cell array containing all created segments of the path.
    % - init_path: initial reference path.
    % - enders: all end points of created path segments.
    % - starters: all start points of created path segments.
    %
    % Outputs:
    % - seg_new: Updated with newly created path segments.
    % - init_path: points used as start and end of segments removed.
    % - enders: all end points of created path segments.
    % - starters: all start points of created path segments.

    seg_add = cell(length(optns)-1,1);
    keep = zeros(length(optns)-1,1);
    parfor ii = 1 : length(optns)-1
        
        % check if waypoint ii can get to point ii+j
        mid_start(ii,:) = optns(ii,:);
        if max(sum(ismember(starters,mid_start(ii,:)),2)) ~= 2
            success = 0;
            j=0;
            ii_limit = 0;
            M=[];
            while success == 0 && ii>=ii_limit && ii+j<size(optns,1)  % keep within segments
                j = j+1;
    
                mid_finish(ii,:) = optns(ii+j,1:2);
                
                % Use the solve_maze_iter function to find the path 
                [path,success,M,~] = solve_maze_iter(Temp,mid_start(ii,:),mid_finish(ii,:),M,2);
    
                if success == 1
                    % save path between optns waypoints by referring to init_path
                    seg_add(ii) = {path};
                    
                    % Iend for optns
                    ii_limit = find(ismember(optns(:,1:2),mid_finish(ii,:),'rows'));
                    
                    % save seg details
                    keep(ii) = 1;
    
                    break
                end
            end
        end
    end    

    % If a new segment is found, update the relevant variables
    if max(keep)==1
        for i = 1 : length(seg_add)
            if keep(i) == 1
                seg_new = [seg_new,seg_add(i)];
                % save previous start and end points
                starters = [starters;mid_start(i,:)];
                enders = [enders;mid_finish(i,:)];

                % remove section with new path
                init_path = rmv_frm_path(init_path,mid_start(i,:),starters,mid_finish(i,:),enders);
            end
        end
    end  
end

function [new_path] = build_newpath(init_path,range,img,ref)
    % Function takes an initial path, a threshold range, an image, 
    % and a reference set of points as inputs, and returns a new 
    % refined path as output. The function goes back through the initial 
    % path and uses image processing to refine the path by removing loops 
    % and connecting disjointed segments.
    % 
    % Inputs:
    % init_path: initial path.
    % range: threshold range for image processing.
    % img: input image.
    % ref: reference points for start and end of path.
    %
    % Outputs:
    % new_path: newly defined tidal channel path.

    diff=rescale(img,1,255);    % Convert image to grayscale and rescale to [1, 255]

    % Initialize variables
    seg_new = {};
    starters = init_path(end,:);
    enders = init_path(1,:);    % first path entry
    mn_sms = 0;
    k_sz = 5;
    thresh = .7*((k_sz*2)^2);


    while ~isempty(init_path)   % While there are still points in the initial path

        % Check which of init_path stay within threshold
        [Temp] = image_varthresh_range(diff, range);
        optns = reachable_points(Temp,init_path);

        % If the mean sum of pixels surrounding optns points is below threshold
        if mn_sms<thresh && ~isempty(optns)
            sms = zeros(size(optns,1),1);
            for z = 1 : size(optns,1)
                [~,k_idx] = create_kernel(k_sz,optns(z,:),Temp);
                sms(z)=sum(Temp(k_idx));
            end
            mn_sms = mean(sms);
        end
        
        % if the mean sum is greater than threshold
        if mn_sms>=thresh
            [seg_new,init_path,enders,starters] = connect_path_par(Temp,optns,seg_new,init_path,enders,starters);
        end
        
        % Increase the range until reaching a maximum and resetting it
        if range > 50
            range=range+25;
        elseif range > 25
            range=range+5;
        else
            range=range+1;
        end

        if size(init_path,1)==1
            init_path=[];
        elseif range>255
            if thresh == 0
                init_path=[];
            end
            range = 1;
            thresh=0;
        end
        if range > 150
            range = 255;
        end

    end

    % Orders path such that the start and end segments are currectly placed
    path=[];
    for i = 1 : length(seg_new)
        path = [path;seg_new{i}];
    end
    dist = sqrt((repmat(ref(:,1)',size(path,1),1)-path(:,1)).^2+(repmat(ref(:,2)',size(path,1),1)-path(:,2)).^2);
    [~,I1]=min(dist(:,1));
    [~,I2]=min(dist(:,2));
    path = [path(I1,:);path;path(I2,:)];

    new_path = remove_loops(path);      % find single path with BFS
end

function [kernel,k_idx] = create_kernel(k_sz,optn,Temp)
    % Function creates a kernel with a square shape of size k_sz x k_sz centred on optn
    %
    % Input:
    %   k_sz: kernel size (in pixels)
    %   optn: (x,y) coordinates of the kernel centre
    %   Temp: input image
    %
    % Output:
    %   kernel: 2D matrix with the coordinates of the kernel pixels
    %   k_idx: vector with the linear indices of the kernel pixels in the input image
    
    x_krnl = optn(:,1)-k_sz:1:optn(:,1)+k_sz;
    y_krnl = optn(:,2)-k_sz:1:optn(:,2)+k_sz;
    
    x_rep = repmat(x_krnl,1,length(y_krnl));
    y_rep= reshape(repmat(y_krnl,length(x_krnl),1),[],1);

    kernel = [x_rep',y_rep];

    k_idx = sub2ind(size(Temp), kernel(:,1),kernel(:,2));
end

function [spaced_path] = downsampling_points(new_path,scale)
    % Function to downsample path points by keeping an evenly space subset
    % of points.
    
    tpath=unique(new_path(:,1:2),'stable','rows');
    dist = ((tpath(1:end-1,1)-tpath(2:end,1)).^2+(tpath(1:end-1,2)-tpath(2:end,2)).^2).^(0.5);
    spaced_path = tpath(1,:);

    % go through tpath removing close proximity points
    combine_dist = 0;
    for i = 2 : length(tpath)
        mn_path = scale*mean(dist);
        if (combine_dist+dist(i-1))>mn_path
            spaced_path=[spaced_path;tpath(i,:)];
            combine_dist = 0;
        else
            combine_dist = combine_dist+dist(i-1);
        end
    end

    if spaced_path(1,1)<spaced_path(end,1)
        spaced_path = flip(spaced_path);
    end
end

function [save_img,file] = create_figure(nm,save_img,foldernm)
    % Function to open image

    file = strcat(foldernm,"VV_only/",nm,"_VV.tif");
    if isfile(file)
        clear geotiffinfo
        [save_img,~,~] = input_process_image(file);
    end
end

function [] = display_path(i,img,old_path,nm)
    % Function to create image, add path overlay, and then save

    img(img==0)=255;
    figure;
    h=imshow(rescale(img,0,1));

    if i>0
        hold on
        scatter(old_path(:,2),old_path(:,1),0.5,'c.')
    end   
 
    saveas(gcf,['Figures/',nm,'.png'])
    close gcf
end

function [] = export_coordinates(spaced_path,lat,lon,nm)
    % Function to export lat lon coordinates to csv

    t=0;
    for ii = 1 : length(spaced_path)
        if ~isempty(spaced_path(ii,2))
            t=t+1;
            path_coord(t,:)=[lat(spaced_path(ii,1),spaced_path(ii,2)),lon(spaced_path(ii,1),spaced_path(ii,2))];
        end
    end
    coordinates = [path_coord(1 : 1 : end,1),path_coord(1 : 1 : end,2)];
    writematrix([[1:length(coordinates)]',coordinates],strcat('Coordinates/Glencaple_',nm(1:10),'.csv'))
end

function [smth_img] = smooth_image(img,windowSize)
    % Function to smooth image using defined windowSize

    kernel = ones(windowSize);
    kernel = kernel / sum(kernel(:));
    smth_img = conv2(img, kernel, 'same');
end

function new_path = river_path_finder(path,img,river_ref)
    % Function applies gradient descent of pixels towards defined river_ref
    new_path=[];
    img(img<0.001)=max(img(:));
    d = [0 0; 0 1; 0 -1; 1 0; -1 0; 1 1; 1 -1;-1 1;-1 -1];
    for i = 1 : size(path,1)
        loc = path(i,:);
        I=0; loop=0;
        while I~=1                                          % while the starting pixel is not selected
            loop=loop+1;
            srch = loc+d;
            srch(srch(:,1)>size(img,1),1)=size(img,1);      %  if srch is outside the bounds of the image Temp
            srch(srch(:,2)>size(img,2),2)=size(img,2);      %  if srch is outside the bounds of Temp
            idx = sub2ind(size(img), srch(:,1), srch(:,2));
            [~,I]=min(abs(img(idx)-river_ref));             % minimum difference from river ref
            loc = srch(I,:);
            if loop > 100                                   % number of pixels transitioned
                break
            end
        end
        new_path=[new_path;loc];
    end

end

function [ref] = image_threshrange(save_img, thresh, range)
    %Function finds navigable pixels based on single intensity value 
    % for ideal conditions image

    ref = zeros(size(save_img,1),size(save_img,2));
    ref(save_img<=thresh+range & save_img>=thresh-range)=1;
end

function [ref] = image_varthresh_range(diff, range)
    % Function finds navigable pixels based on manipulated image

    ref = zeros(size(diff,1),size(diff,2));
    ref(abs(diff)<=range)=1;
end

function [PDImg] = river_ref_diff(save_img,riv_rollmn)
    % Function calculates the difference between the input image
    % and a reference image based on the pixel proximity to path points and
    % the path point intensities.
    % 
    % Inputs:
    % - save_img: input image.
    % - riv_rollmn: rolling intensity average associated with channel path points.
    %
    % Outputs:
    % - diff: Path difference image
    
    % create coordinates matrix to calculate reference image
    A=repmat(1:size(save_img,1),size(save_img,2),1)';
    B=repmat(1:size(save_img,2),size(save_img,1),1);
    coords = [A(:),B(:)];

    % calculate the reference image
    [var_ref,~] = var_river_ref_repeat(coords,riv_rollmn);
    ref_mat = reshape(var_ref,size(save_img,1),size(save_img,2));

    PDImg = abs(save_img-ref_mat);  % path difference image
end

function [var_ref,dist_ref] = var_river_ref_repeat(old_path,riv_rollmn)
    % Function to find distance from riv_rollmn points to all other image pixels

   dist = zeros(size(old_path,1),size(riv_rollmn,1));
    for i = 1 : size(riv_rollmn,1)
        dist(:,i) = sqrt((old_path(:,1)'-riv_rollmn(i,1)).^2+(old_path(:,2)'-riv_rollmn(i,2)).^2);
    end
    [~,I_sav] = min(dist,[],2);

    var_ref = riv_rollmn(I_sav,3);

    [~,indA]=min(dist,[],1);
    indB=sub2ind(size(dist),indA,(1:size(dist,2)));
    dist_ref=dist(indB);
end

function [mean_channel] = mean_channel_intensity(spaced_path,img)
    % Function to find mean channel intensity for points on spaced_path

    sar_v = [];
    d = [0 0; 0 1; 0 -1; 1 0; -1 0; 1 1; 1 -1;-1 1;-1 -1];
    
    for i = 1 : length(spaced_path)
        pix = [spaced_path(i,1),spaced_path(i,2)]+d;

        pix(pix(:,1)>size(img,1),1)=size(img,1);    % ensure window is within bounds
        pix(pix(:,2)>size(img,2),2)=size(img,2);

        mn_val = [];
        for j = 1 : length(pix)
            mn_val = [mn_val,img(pix(j,1),pix(j,2))];
        end
        sar_v = [sar_v;mn_val];
    end
    mean_channel = mean(sar_v(~isnan(sar_v)));

end

function [check,centre] = check_SAR(img,riv_rollmn)
    % Function to check if riv_rollmn points (corresponding to the tidal 
    % channel) are the lowest intensity feature (associated with the 
    % lowest intensity peak in the image)

    I = sub2ind(size(img),riv_rollmn(:,1),riv_rollmn(:,2));
    vals = img(I);
    mn_val = median(vals(~isnan(vals)));

    figure;
    h=histogram(img(~isnan(img)));
    h.BinWidth = 5;
    title(num2str(mn_val))

    % Retrieve some properties from the histogram
    V = h.Values;
    E = h.BinEdges;

    peak = find(islocalmax(V,'MinSeparation',3,'MinProminence',500));
  
    centre = E(peak(1))+2.5;

    if centre+10 < mn_val
        check = 1;
    else
        check = 0;
    end
    close gcf
end

function [] = saving_path(nm,file,new_path,mean_channel,cat_all,first_path,folder_sv)
    % Function that saves key variables

    sav_dat.name = file;
    sav_dat.new_path = new_path;
    sav_dat.first_path = first_path;
    sav_dat.mean_channel = mean_channel;

    if size(cat_all,3)<5    % If the number of previous masks is less than 5, average all the masks
        sav_dat.cat_all = cat_all(:,:,:); 
    else
        sav_dat.cat_all = cat_all(:,:,end-4:end); 
    end
    
    save(strcat(folder_sv,nm,'.mat'),'-struct','sav_dat', '-v7.3')
end

function [spaced_path] = fix_spaced_path(i,spaced_path,old_path,x_def,y_def,takePath,first)
    % Function that extends spaced_path with previous path, when the image
    % for spaced_path does not cover the whole tidal channel

    [spaced_path] = check_for_flip(spaced_path,y_def,x_def);
    
    if takePath > 0 && i>first
        dist_old_1 = dist_compare([old_path(1,1),old_path(1,2)],[y_def(2),x_def(2)]);
        dist_old_end = dist_compare([old_path(end,1),old_path(end,2)],[y_def(2),x_def(2)]);
        if dist_old_1 > dist_old_end                                        % if old_path end is nearer to route end then flip
            old_path = flip(old_path);
        end
        if takePath > 1
            spaced_path = [spaced_path;flip(old_path(1:takePath,:))];       % take the old_path nearest route end
        end
    end

    spaced_path = [y_def(1),x_def(1);spaced_path;y_def(end),x_def(end)];    % Add ref start and end
end

function [path] = check_for_flip(path,y_def,x_def)
    % Function to check if path order needs reversed to be in the format 
    % expected by the algorithm

    dist_sp_1 = dist_compare([path(1,1),path(1,2)],[y_def(2),x_def(2)]);
    dist_sp_end = dist_compare([path(end,1),path(end,2)],[y_def(2),x_def(2)]);

    if dist_sp_1 < dist_sp_end  % if spaced path start is nearer to route end then flip
        path = flip(path);
    end
end