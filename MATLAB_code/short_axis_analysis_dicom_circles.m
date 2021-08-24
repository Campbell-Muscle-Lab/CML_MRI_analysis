function short_axis_analysis

% Variables

dicom_file_string = '../data/enIm5.dcm';
frame_number=8;

% Gauss filter SD
gauss_filter_sd = 1;

% Radius range for initial search
radius_range = [10 100];
n_circles = 1;
circle_expansion = 2;   % big enough to include both ventricles

% Adaptive contouring
contour_contrast = [0 0.5]; % Contour contrast enhancement
no_of_iterations = 100;
smooth_factor = 0.5;
contraction_bias = -0.5;

% Display
zoom_limits = [50 150];

% Angles and thresholds to scan for septal and opposite epicardial walls
septal_angles = linspace(-60, 20, 10);
peak_prominence = 0.1;

opposite_angles = linspace(140, 180, 10);
opposite_threshold = 0.4;

interp_x = 0.1;

% Code

% Read dicom
dic = dicomread(dicom_file_string);
[x_pixels, y_pixels, ~, no_of_frames] = size(dic);

% Create a figure
analysis_figure = figure(1);
clf;
no_of_rows = 4;
no_of_cols = 4;

% Loop through frames
for frame_counter = frame_number : frame_number
    
    % Load frame as double and normalize
    im_raw = double(dic(:,:,frame_counter));
    im_raw = im_raw./max(im_raw(:));
    
    % Smooth
    if (gauss_filter_sd > 0)
        im_f = imgaussfilt(im_raw, gauss_filter_sd);
    else
        im_f = im_raw;
    end
    
    % Find circles
    [centers, radii] = imfindcircles(im_f, radius_range,'ObjectPolarity','bright');
    
    % Get mask that encompasses both ventricles
    im_circle = zeros(size(im_f));
    for i = 1:x_pixels
        for j = 1:x_pixels
            h = hypot(i-centers(1,1),j-centers(1,2));
            if (h<(circle_expansion*radii(1)))
                im_circle(j,i)=1;
            end
        end
    end
    
    % Calculate std of image, adaptive binarize, and thicken edges
    % This should enclose LV and RV in white rings
    im_std = imadjust(stdfilt(im_f,true(3)));
    im_std_2 = imbinarize(im_std,'adaptive');
    im_std_3 = imclose(im_std_2,strel('disk',1));
    
    % Restrict edge picture to circular mask
    im_std_3(~im_circle)=0;
    % Fill in from edges so that we just have the LV and RV rings
    im_filled_ventricles = imcomplement(imfill(im_std_3,[1 1]));
    % Keep 2 largest components, should be RV (top) and LV (bottom)
    im_filled_ventricles = bwareafilt(im_filled_ventricles,2);
    im_filled_ventricles = imfill(im_filled_ventricles,'holes');
    
    % Find LV and RV
    im_L = bwlabel(im_filled_ventricles);
    % Adapt image for active contouring
    im_cont = imcomplement(im_f);
    im_cont = imadjust(im_cont, contour_contrast, [0 1]);
    im_cont(~im_circle)=1;
    
    % Now do active contouring for both RV (1) and LV (2)
    for i = 1 : 2
        im_m = zeros(size(im_L));
        im_m(im_L==i)=1;
        im_mask{i} = im_m;
        im_c{i} = activecontour(im_cont, im_mask{i}, no_of_iterations, ...
            'Chan-Vese', ...
            'SmoothFactor', smooth_factor, ...
            'ContractionBias', contraction_bias);
    end
    
    % Find the centroid
    lv_props = regionprops(im_c{2}, {'Centroid'});
    temp = cat(1, lv_props);
    lv_centroid = round(temp(1).Centroid);

    % Display frame images
    subplot_counter = display_frame_images()

    % Loop through septal angles
    for septal_counter = 1:numel(septal_angles)
        
        % Rotate the image around the centroid
        im_rot = rotateAround(im_f, lv_centroid(2), lv_centroid(1), ...
                    septal_angles(septal_counter));
                
        [sr(septal_counter), out] = profile_analysis('septum');
        
        display_angle_images(subplot_counter + 1, 'septum', out);
    end
    
    % Loop through opposite angles
    for opposite_counter = 1:numel(opposite_angles)
        
        % Rotate the image around the centroid
        im_rot = rotateAround(im_f, lv_centroid(2), lv_centroid(1), ...
                    opposite_angles(opposite_counter));
                
        [or(opposite_counter), out] = profile_analysis('opposite');
        
        display_angle_images(subplot_counter + 3, 'opposite', out);
    end
    
    figure(2);
    clf
    subplot(1,3,1);
    imagesc(im_raw);
    title('Raw');
   
    subplot(1,3,2);
    imagesc(im_f);
    title(sprintf('Filtered with Gaussian SD of %.2f',gauss_filter_sd));
    
    ax = subplot(1,3,3)
    colormap(gray);
    hold on;
    g = repmat(zeros(size(im_f)),[1 1 3]);
    g(:,:,2)=1;
    imagesc(im_f);
    title('Analysis result');
    h = image(g);
    set(h, 'AlphaData', 0.5*(im_c{2} + im_c{1}));
    plot(lv_centroid(1), lv_centroid(2), 'ro');
    for septal_counter = 1 : numel(septal_angles)
        sx(septal_counter) = lv_centroid(1) + sr(septal_counter)*sind(septal_angles(septal_counter))
        sy(septal_counter) = lv_centroid(2) - sr(septal_counter)*cosd(septal_angles(septal_counter))
    end
    plot(sx, sy, 'bo');

    for opposite_counter = 1 : numel(opposite_angles)
        ox(opposite_counter) = lv_centroid(1) + or(opposite_counter)*sind(opposite_angles(opposite_counter))
        oy(opposite_counter) = lv_centroid(2) - or(opposite_counter)*cosd(opposite_angles(opposite_counter))
    end
    plot(ox, oy, 'go');
    
    % Coordinates for ellipse
    xe = [sx ox];
    ye = [sy oy];
    
    e_out = fit_ellipse(xe',ye')
    
    % Draw the ellipse
    R = [ e_out.cos_phi e_out.sin_phi ; -e_out.sin_phi e_out.cos_phi ];
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = e_out.X0 + e_out.a*cos(theta_r);
    ellipse_y_r     = e_out.Y0 + e_out.b*sin(theta_r);
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    plot(rotated_ellipse(1,:), rotated_ellipse(2,:), 'y-', 'LineWidth',2);
  
    xlim([1 x_pixels]);
    ylim([1 y_pixels]);
    set(gca, 'YDir', 'reverse');


    
    if (1)
        break
    end
    
    if (video_output)
        close(v);
    end
end
        % Nested functions
        function [ri, out] = profile_analysis(region_string)
            % Returns epicardial radius

            % Pull off line profile
            p_r = 1:lv_centroid(2);
            p_z = im_rot(p_r, lv_centroid(1));
            p_z(end) = 0;
            
            % Interpolate
            p_ri = 1:interp_x:p_r(end);
            p_zi = interp1(p_r, p_z, p_ri, 'spline');
            
            % Find peaks and troughs
            [~, p_loci]=findpeaks(p_zi, 'MinPeakProminence', ...
                            peak_prominence*max(p_zi));
            [~, t_loci]=findpeaks(-p_zi, 'MinPeakProminence', ...
                            peak_prominence*max(p_zi));
            
            switch region_string
                case 'septum'
                    th = mean(p_zi([p_loci(end) t_loci(end)]));
                    indi = find((p_zi>th)&(p_ri < p_ri(t_loci(end))), 1, 'last');
                    ri = p_ri(numel(p_zi)-indi);
                case 'opposite'
                    % Set a threshold for the far side of the LV wall
                    th = opposite_threshold * p_zi(p_loci(end));
                    indi = find((p_zi<th)&(p_ri < p_ri(p_loci(end))), 1, 'last');
                    ri = p_ri(numel(p_zi)-indi);
            end
            
            out.p_r = p_r;
            out.p_z = p_z;
            out.p_ri = p_ri;
            out.p_zi = p_zi;
            out.p_loci = p_loci;
            out.t_loci = t_loci;
            out.th = th;
            out.indi = indi;
            
        end

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
            title('std - Standard deviation');
            colorbar;

            % Frangi ventricles
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_std_2);
            title('std2 - Opened to thicken');
            colorbar;


            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_std_3);
            title('std3 - Restricted to circle');
            colorbar;
            
            % ventricles
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_filled_ventricles);
            colorbar;
            title('LV and RV seeds')

            
            % Edges
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            imagesc(im_L);
            colorbar;
            title('LV and RV labels');

            % Edges
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            imagesc(im_cont);
            colorbar;
            title('Image for adaptive contouring');

            
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
            xlim([zoom_limits(1) zoom_limits(2)]);
            ylim([zoom_limits(1) zoom_limits(2)]);
            colorbar;
            title('Ventricles');
            
        end

        function display_angle_images(subplot_counter, sub_title, li)

            % Rotated frame
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            imagesc(im_rot);
            colorbar;
            hold on;
            plot(lv_centroid(1), lv_centroid(2), 'ro');
            plot(lv_centroid(1)*[1 1], [1 lv_centroid(2)], 'b-');
            title([sub_title '-rotated frame']);
 
            % Plot the line profile
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            colormap(gray);
            cla;
            hold on;
            plot(li.p_r, li.p_z, 'b-');
            plot(li.p_ri, li.p_zi, 'm-');
            plot(li.p_ri(li.p_loci), li.p_zi(li.p_loci), 'go');
            plot(li.p_ri(li.t_loci), li.p_zi(li.t_loci), 'ro');
            plot(li.p_r([1 end]), li.th*[1 1], 'c-');
            plot(li.p_ri(li.indi), li.p_zi(li.indi), 'md');
            title([sub_title '-line profile']);
        end
end
