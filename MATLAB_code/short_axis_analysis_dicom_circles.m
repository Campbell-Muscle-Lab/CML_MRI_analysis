function short_axis_analysis

% Variables

dicom_file_string = '../data/enIm10.dcm';
frame_number=8;

radius_range = [10 100];
n_circles = 1;

ll = [50 150];

lv_angles = linspace(-60, 20, 10);

rv_angles = linspace(140, 200, 10);

% Code

% Read dicom
dic = dicomread(dicom_file_string);
[x_pixels, y_pixels, ~, no_of_frames] = size(dic);

% Create a figure
analysis_figure = figure(1);
clf;
no_of_rows = 4;
no_of_cols = 4;

wall_thickness = NaN*ones(no_of_frames, 1);

% Loop through frames
for frame_counter = frame_number : frame_number
    
    % Load frame as double and normalize
    im_f = double(dic(:,:,frame_counter));
    im_f = im_f./max(im_f(:));
    
    % Find circles
    [centers, radii] = imfindcircles(im_f, radius_range,'ObjectPolarity','bright');
    
    % Get mask that encompasses both ventricles
    im_circle = zeros(size(im_f));
    for i = 1:x_pixels
        for j = 1:x_pixels
            h = hypot(i-centers(1,1),j-centers(1,2));
            if (h<(3*radii(1)))
                im_circle(j,i)=1;
            end
        end
    end
    
    im_std = imadjust(stdfilt(im_f,true(3)));
    im_std_2 = imbinarize(im_std,'adaptive');
    im_std_3 = imclose(im_std_2,strel('disk',1));
    
    % Image ventricles
    im_ventricles = im_std_3;
    im_ventricles(~im_circle)=0;
    
    im_filled_ventricles = imcomplement(imfill(im_ventricles,[1 1]));
    im_filled_ventricles = bwareafilt(im_filled_ventricles,2);
    
    % Find LV and RV
    im_L = bwlabel(im_filled_ventricles);
    
    im_cont = imcomplement(imadjust(im_f,[0.5 1],[0 1]));
%     n=3;
%     se = ones(n); se = se/sum(se(:));
%     im_cont = conv2(im_cont, se, 'same');
    
    ni=100;
    for i = 1 : 2
        im_m = zeros(size(im_L));
        im_m(im_L==i)=1;
        im_mask{i} = im_m;
        im_c{i} = activecontour(im_cont, im_mask{i}, ni, 'Chan-Vese', ...
            'SmoothFactor',0.5, ...
            'ContractionBias',-0.5);
    end
    
    lv_props = regionprops(im_c{2}, {'Centroid'});
    temp = cat(1, lv_props)
    lv_centroid = round(temp(1).Centroid)

    % Display frame images
    subplot_counter = display_frame_images()

    % Loop through lv angles
    for lv_counter = 1:numel(lv_angles)
        im_rot = rotateAround(im_f, lv_centroid(2), lv_centroid(1), lv_angles(lv_counter));
        
        p_r = 1:lv_centroid(2);
        p_z = im_rot(p_r, lv_centroid(1));
        p_z(end)=0;
        
        [~, p_loc]=findpeaks(p_z, 'MinPeakProminence',0.1*max(p_z));
        p_loc = p_loc(end-1:end);
        [~, t_loc]=findpeaks(-p_z, 'MinPeakProminence',0.1*max(p_z));
        t_loc = t_loc(end);
        
        th = mean(p_z([p_loc(end) t_loc(end)]));
        
        ind = find((p_z>th)&(p_r < t_loc(end))', 1, 'last')
        lz(lv_counter) = numel(p_z)-ind
       
        display_angle_images(subplot_counter+1);
    end
    
    % Loop through rv angles
    for rv_counter = 1:numel(rv_angles)
        im_rot = rotateAround(im_f, lv_centroid(2), lv_centroid(1), rv_angles(rv_counter));
        
        p_r = 1:lv_centroid(2);
        p_z = im_rot(p_r, lv_centroid(1));
        p_z(end)=0;
        
        [~, p_loc]=findpeaks(p_z, 'MinPeakProminence',0.1*max(p_z));
        th = 0.25 * p_z(p_loc(end));
        
        ind = find((p_z<th)&(p_r < t_loc(end))', 1, 'last')
        rz(lv_counter) = numel(p_z)-ind
       
%         display_angle_images(subplot_counter+1);
    end
    
    figure(2);
    clf
    subplot(1,2,1);
    imagesc(im_f);
   
    ax = subplot(1,2,2)
    colormap(gray);
    hold on;
    g = repmat(zeros(size(im_f)),[1 1 3]);
    g(:,:,2)=1;
    imagesc(im_f);
    h = image(g);
    set(h, 'AlphaData', 0.5*(im_c{2} + im_c{1}));
    plot(lv_centroid(1), lv_centroid(2), 'ro');
    for lv_counter = 1 : numel(lv_angles)
        lx(lv_counter) = lv_centroid(1) + lz(lv_counter)*sind(lv_angles(lv_counter))
        ly(lv_counter) = lv_centroid(2) - lz(lv_counter)*cosd(lv_angles(lv_counter))
    end
    plot(lx, ly, 'bo');

    for rv_counter = 1 : numel(rv_angles)
        rx(rv_counter) = lv_centroid(1) + rz(lv_counter)*sind(rv_angles(rv_counter))
        ry(rv_counter) = lv_centroid(2) - rz(lv_counter)*cosd(rv_angles(rv_counter))
    end
    plot(rx, ry, 'go');
    
    xe = [lx rx];
    ye = [ly ry];
    
    fit_ellipse(xe',ye',ax)
    
    xlim([1 x_pixels]);
    set(gca, 'YDir', 'reverse');


    
    if (1)
        break
    end
    
    if (video_output)
        close(v);
    end
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

            % Circles
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_f);
            hold on;
            viscircles(centers(1:n_circles,:), radii(1:n_circles), 'EdgeColor','b');
            title('Circles');
            colorbar;
            
            % Mask
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(imoverlay(im_f, im_circle, 'green'));
            title('mask');
            colorbar;

            % Ventricles
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_std);
            title('std');
            colorbar;

            % Frangi ventricles
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_std_2);
            title('std2');
            colorbar;


            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_std_3);
            title('std3');
            colorbar;

            % ventricles
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_ventricles);
            colorbar;
            
            % ventricles
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_filled_ventricles);
            colorbar;

            
            % Edges
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            imagesc(im_L);
            colorbar;
            
            % Edges
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            imagesc(im_cont);
            hold on;
            cm = jet(3);
            for i=1:2
                visboundaries(im_c{i},'Color',cm(i,:));
            end
            plot(lv_centroid(1), lv_centroid(2), 'ro');
            xlim([ll(1) ll(2)]);
            ylim([ll(1) ll(2)]);
            colorbar;
            title('Ventricles');
            
        end

        function display_angle_images(subplot_counter)

            % Rotated frame
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_rot);
            colorbar;
            hold on;
            plot(lv_centroid(1), lv_centroid(2), 'ro');
%             plot(lv_centroid(1)*[1 1], [1 lv_centroid(2)], 'r-');
%             plot(lv_centroid(1), x_out, 'bo');
%             plot(lv_centroid(1), x_in, 'go');
            title('Rotated frame');
% 
            % Plot the line profile
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            hold on;
            plot(p_r, p_z, 'b-');
            plot(p_r(p_loc), p_z(p_loc), 'go');
            plot(p_r(t_loc), p_z(t_loc), 'ro');
            plot(p_r([1 end]), th*[1 1], 'c-');
            plot(p_r(ind), p_z(ind), 'mo');
            
%             hold on;
%             plot(p_x([1 end]), wall_thresh*[1 1], 'c-');
%             plot(x_last_peak, p_y(x_last_peak), 'co');
%             plot(x_last_trough, p_y(x_last_trough), 'co');
%             plot(x_out, p_y(x_out), 'bo');
%             plot(x_in, p_y(x_in), 'go');
            title('Line profile');
%             
%             % Plot a zoom
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             colormap(gray);
%             cla;
%             imagesc(im_rot);
%             colorbar;
%             hold on;
%             plot(lv_centroid(1), lv_centroid(2), 'r+');
%             plot(lv_centroid(1)*[1 1], [1 lv_centroid(2)], 'r-');
%             plot(lv_centroid(1), x_out, 'bo');
%             plot(lv_centroid(1), x_in, 'go');
%             x_limits = lv_centroid(1) + 2 * (lv_centroid(1) - x_out) * [-1 1];
%             x_limits(x_limits<1)=1;
%             x_limits(x_limits>x_pixels)=x_pixels;
%             xlim(x_limits);
%             y_limits = lv_centroid(2) + 2 * (lv_centroid(2) - x_out) * [-1 1];
%             y_limits(y_limits<1)=1;
%             y_limits(y_limits>y_pixels)=y_pixels;
%             ylim(y_limits);
%             title('Zoom');
%             
%             % Plot a zoom
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             hold on;
%             cm = paruly(no_of_frames);
%             plot(ang, wall_thickness(frame_counter, :), '-o', ...
%                 'Color', cm(frame_counter, :));
%             xlabel('Angle');
%             ylabel('Wall thickness');            
        end
end
