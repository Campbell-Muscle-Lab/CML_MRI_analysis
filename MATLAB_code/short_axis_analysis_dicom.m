function short_axis_analysis

% Variables
dicom_file_string = 'c:\temp\dicom.dcm';
no_of_segmentation_levels = 6;
no_of_angles = 50;
min_peak_depth = 0.025;

wall_prop_thresh = 0.5;
 
v = VideoWriter('example_analysis_2', 'Motion JPEG 2000');
v.FrameRate = 3;
open(v);

% Code

% Read dicom
dic = dicomread(dicom_file_string);
[x_pixels, y_pixels, ~, no_of_frames] = size(dic);

% Create a figure
analysis_figure = figure(1);
clf;
no_of_rows = 3;
no_of_cols = 4;

wall_thickness = NaN*ones(no_of_frames, no_of_angles);

% Loop through frames
for frame_counter = 11 : no_of_frames
    
    % Load frame as double
    im_f = double(dic(:,:,frame_counter));
    
    % Enhance contrast
    im_f_enhanced = (im_f./max(im_f(:)));
        
    % Threshold the image and find the brightest level
    im_levels = imquantize(im_f_enhanced, ...
        multithresh(im_f_enhanced, no_of_segmentation_levels));
    im_brightest = zeros(size(im_levels));
    im_brightest(im_levels==(no_of_segmentation_levels+1)) = 1;
    
    % Watershed
    im_distance = -bwdist(~im_brightest);
    im_L_distance = watershed(im_distance);
    im_L_distance(~im_brightest)=0;
    
    % Find LV as the biggest region
    s = regionprops(bwlabel(im_L_distance), {'centroid', 'PixelIDxList', 'area'})
    
    s_areas = cat(1, s.Area);
    [~, si] = sort(s_areas, 'descend');
    lv_pixels_idx = s(si(1)).PixelIdxList;
    lv_centroid = round(s(si(1)).Centroid);
    
    % Set the lv_seed
    im_lv_seed = zeros(size(im_f));
    im_lv_seed(lv_pixels_idx)=1;
    
    % Fill the LV
    im_lv_filled = im_f_enhanced;
    im_lv_filled(lv_pixels_idx) = mean(im_f_enhanced(lv_pixels_idx));

    % Display frame images
    subplot_counter = display_frame_images()
    
    % Scan through rotations
    ang = linspace(0, 360, no_of_angles);
    for angle_counter = 1 : no_of_angles
        % Rotate the image
        im_rot = rotateAround(im_lv_filled, ...
            lv_centroid(2), lv_centroid(1), ang(angle_counter));
      
        % Get a line profile
        c = 1:lv_centroid(2);
        ro = lv_centroid(1);
        p_y = im_rot(c, ro);
        p_x = 1:numel(p_y);
        
        % Find peaks and troughs
        pd = find_peaks('x', p_x, 'y', p_y, ...
                'min_rel_delta_y', min_peak_depth, ...
                'min_x_index_spacing', 1);

        % Last trough
        x_last_trough = pd.min_indices(end);
        % Last peak
        x_last_peak = pd.max_indices(end);
        % Thresh is wall_threshold betwee between
        wall_thresh = double(p_y(x_last_trough) + ...
            wall_prop_thresh * (p_y(x_last_peak) - p_y(x_last_trough)))
        
        % Find last point
        x_out = find((p_y > wall_thresh) & (p_x < x_last_trough)', ...
                    1, 'last');
        x_in = find((p_y > wall_thresh) & (p_x > x_last_trough)', ...
                    1, 'first');

        % Deduce the wall thickness and save in array
        wall_thickness(frame_counter, angle_counter) = [x_in - x_out];                

        display_angle_images(subplot_counter+1);
        drawnow;
        pause(0.2);
        
        fr = getframe(gcf);
        writeVideo(v,fr);

    end
    
    if (1)
        break
    end
    
    close(v);
end

        % Nested function
        function subplot_counter = display_frame_images()

            % Display
            figure(analysis_figure)
            clf
            subplot_counter = 1;

            % Raw frame
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_f);
            colorbar;
            title('Raw frame');

            % Contrast enchanced
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_f_enhanced);
            colorbar;
            title('Contrast enhanced');

            % Segmentation
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_levels);
            title('Threshold');
            colorbar;

            % Brightest segment
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_brightest);
            title('Brightest segment');

            % Watershed
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_distance);
            title('Distance image');

            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_L_distance);
            title('Watershed');

            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_lv_seed);
            title('LV segment');

            % Fill the LV
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_lv_filled);
            hold on;
            plot(lv_centroid(1), lv_centroid(2), 'r+');
            title('LV filled');
        end

        function display_angle_images(subplot_counter)

            % Rotated frame
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_rot);
            colorbar;
            hold on;
            plot(lv_centroid(1), lv_centroid(2), 'r+');
            plot(lv_centroid(1)*[1 1], [1 lv_centroid(2)], 'r-');
            plot(lv_centroid(1), x_out, 'bo');
            plot(lv_centroid(1), x_in, 'go');
            title('Rotated frame');

            % Plot the line profile
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            plot(p_x, p_y, 'r-');
            hold on;
            plot(p_x([1 end]), wall_thresh*[1 1], 'c-');
            plot(x_last_peak, p_y(x_last_peak), 'co');
            plot(x_last_trough, p_y(x_last_trough), 'co');
            plot(x_out, p_y(x_out), 'bo');
            plot(x_in, p_y(x_in), 'go');
            title('Line profile');
            
            % Plot a zoom
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_rot);
            colorbar;
            hold on;
            plot(lv_centroid(1), lv_centroid(2), 'r+');
            plot(lv_centroid(1)*[1 1], [1 lv_centroid(2)], 'r-');
            plot(lv_centroid(1), x_out, 'bo');
            plot(lv_centroid(1), x_in, 'go');
            x_limits = lv_centroid(1) + 2 * (lv_centroid(1) - x_out) * [-1 1]
            x_limits(x_limits<1)=1;
            x_limits(x_limits>x_pixels)=x_pixels;
            xlim(x_limits);
            y_limits = lv_centroid(2) + 2 * (lv_centroid(2) - x_out) * [-1 1]
            y_limits(y_limits<1)=1;
            y_limits(y_limits>y_pixels)=y_pixels;
            ylim(y_limits);
            title('Zoom');
            
            % Plot a zoom
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            hold on;
            cm = paruly(no_of_frames);
            plot(ang, wall_thickness(frame_counter, :), '-o', ...
                'Color', cm(frame_counter, :));
            xlabel('Angle');
            ylabel('Wall thickness');            
        end
end
