function short_axis_analysis

% Variables

dicom_file_string = '../data/enIm9.dcm';
frame_number=8;

radius_range = [5 100];
n_circles = 1;

ll = [50 150];

% Code

% Read dicom
dic = dicomread(dicom_file_string);
[x_pixels, y_pixels, ~, no_of_frames] = size(dic);

% Create a figure
analysis_figure = figure(1);
clf;
no_of_rows = 4;
no_of_cols = 4;

% wall_thickness = NaN*ones(no_of_frames, no_of_angles);

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
    
    im_d = imdilate(im_filled_ventricles,strel('disk',3));
    rp = regionprops(im_d,'PixelIdxList');
    im_cv = zeros(size(im_f));
    im_cv(rp(1).PixelIdxList)=1;
    im_cv2 = im_cv;
    im_cv2(im_mask{1}==1)=0;
    im_cv2(im_mask{2}==1)=0;
    
    mm = mean(im_f(im_cv2>0));
    
    
    im_v = im_f;
    im_v(im_mask{1}>0)=mm;
    im_v = imadjust(im_v,[0 1],[0 1]);
% 
%     % Find blob below lv
%     rp = regionprops(im_L, {'centroid'});
%     cent = cat(1,rp.Centroid);
%     lv_cent = round(cent(2,:))
%     
%     im_std_2_circle = im_std_2;
%     im_std_2_circle(~im_circle)=0;
%     im_std_2_circle = bwareafilt(im_std_2_circle,[10 inf]);
%     im_L2 = bwlabel(im_std_2_circle);
%     % Find blob below LV
%     p = im_L2(1:end, lv_cent(1));
%     p2 = p(lv_cent(2):end);
%     p2 = p2(p2>0)
%     p2(diff(p2)==0)=[];
%     
%     im_bot = im_L2;
%     im_bot(im_L2~=p2(2))=0;
%     
%     % Superpose bottom
%     im_bot_2 = imoverlay(im_f, im_bot, 'g');
%     
    
    
%     
%     % Make lv_seed
%     ss = size(im_f)
%     
%     im_mask = zeros(size(im_f));
%     c = round(centers(1,:))
%     im_mask(c(2),c(1))=1;
%     im_mask = imdilate(im_mask,strel('disk',10));
%     
%     ni = 100;
%     im_contour = im_f;
%     im_contour = 100*imadjust(im_contour, [0.0 1],[0 1]);
%     im_ac = activecontour(im_contour, im_mask, ni, 'Chan-Vese', ...
%                 'SmoothFactor',0.5);

    % Display frame images
    subplot_counter = display_frame_images()
    
    % Scan through rotations
%     ang = linspace(0, 360, no_of_angles);
%     for angle_counter = 1 : no_of_angles
%     end
%     
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
            xlim([ll(1) ll(2)]);
            ylim([ll(1) ll(2)]);
            colorbar;
            title('Ventricles');
            
            % CV
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            imagesc(im_cv);
            colorbar;
            
             % CV
            subplot_counter = subplot_counter + 1;
            subplot(no_of_rows, no_of_cols, subplot_counter);
            imagesc(im_v);
            xlim([ll(1) ll(2)]);
            ylim([ll(1) ll(2)]);
            colorbar;
%             
%             % L2 profile
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             plot(p);
%             colorbar;
%             
%             % L2 profile
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             imagesc(im_bot);
%             colorbar;
%   
%              % L2 profile
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             imagesc(im_bot_2);
%             colorbar;
%             
            
            
            
            
            % 
%             % Ventricles
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             colormap(gray);
%             cla;
%             imagesc(im_v)
%             title('Ventricles');
%             colorbar;
% 
%             % Ventricles
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             cla;
%             imagesc(im_l)
%             title('Ventricles');
%             colorbar;
% 
%             % Ventricles
%             subplot_counter = subplot_counter + 1;
%             subplot(no_of_rows, no_of_cols, subplot_counter);
%             cla;
%             imagesc(im_f_contour);
% %             hold on;
% %             cm = jet(3);
% %             for ci=1:(max(im_l(:))+1)
% %                 visboundaries(im_ac{ci}, 'Color', cm(ci,:));
% %             end
%             title('Ventricles');
%             colorbar;
% 
%             
% %             % Contour
% %             subplot_counter = subplot_counter + 1;
% %             subplot(no_of_rows, no_of_cols, subplot_counter);
% %             colormap(gray);
% %             cla;
% %             imagesc(im_f_2);
% %             hold on;
% %             visboundaries(im_c_2, 'Color', 'r');
% %             visboundaries(im_c_3, 'Color', 'b');
% %             visboundaries(im_c_4, 'Color', 'g');
% %             title('Contour');
% %             colorbar;

            
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
            x_limits = lv_centroid(1) + 2 * (lv_centroid(1) - x_out) * [-1 1];
            x_limits(x_limits<1)=1;
            x_limits(x_limits>x_pixels)=x_pixels;
            xlim(x_limits);
            y_limits = lv_centroid(2) + 2 * (lv_centroid(2) - x_out) * [-1 1];
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
