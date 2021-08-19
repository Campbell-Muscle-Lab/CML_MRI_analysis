function short_axis_analysis

% Variables
movie_file_string = '../data/short_axis_mri/lvsa.gif';
lv_seed = [330 240];
no_of_rotations = 10;
wall_prop_thresh = 0.5;

v = VideoWriter('example_analysis_2', 'Motion JPEG 2000');
v.FrameRate = 3;
open(v);

% Code

% Read movie
m = imread(movie_file_string, 'frames','all');
no_of_frames = size(m, 4)

% Create a figure
figure(1);
clf;
no_of_rows = 2;
no_of_cols = 3;

wall_thickness = NaN*ones(no_of_frames, no_of_rotations);

% Loop through frames
for i = 1 : no_of_frames
    i = i
    f = m(:,:,i);
    size(f)
    
    % Threshold the image and fill holes
    im_bw = imbinarize(f);
    im_bw = imfill(im_bw, 'holes');
    
    % Label the image and find the blob with the LV
    im_label = bwlabel(im_bw);
    % Filter the labeled image to show only the LV chamber
    lv_blob = im_label(lv_seed(2), lv_seed(1))
    im_lv_bw = im_bw;
    im_lv_bw(im_label ~= lv_blob) = 0;
    
    % Find the centroid and pixels in the LV
    s = regionprops(im_lv_bw, {'centroid', 'PixelIdxList'});
    lv_centroid = round(cat(1, s.Centroid));
    lv_pixel_ids = cat(1, s.PixelIdxList);
    lv_mean = mean(f(lv_pixel_ids));
    
    % Normalize the LV
    im_lv_filtered = f;
    im_lv_filtered(lv_pixel_ids) = lv_mean;
       
    % Scan through rotations
    ang = linspace(0,360,no_of_rotations)
    for r = 1 : no_of_rotations
        
        im_rot = rotateAround(im_lv_filtered, ...
            lv_centroid(2), lv_centroid(1), ang(r));
        
        % Get a line profile
        c = 1:lv_centroid(2);
        ro = lv_centroid(1);
        p_y = im_rot(c, ro);
        p_x = 1:numel(p_y);
        
        % Find peaks and troughs
        pd = find_peaks('x', p_x, 'y', p_y, ...
                'min_rel_delta_y', 0.05, ...
                'min_x_index_spacing', 5);

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
                
        wall_thickness(i, r) = [x_in - x_out];                
    
        subplot(no_of_rows, no_of_cols, 1);
        colormap(gray);
        cla;
        imagesc(f);
        hold on;
        plot(lv_seed(1), lv_seed(2), 'r+');

        subplot(no_of_rows, no_of_cols, 2);
        cla
        imagesc(im_bw);

        subplot(no_of_rows, no_of_cols, 3);
        cla
        imagesc(im_lv_bw);

        subplot(no_of_rows, no_of_cols, 4);
        cla
        imagesc(im_rot);
        hold on;
        plot(lv_centroid(1), lv_centroid(2), 'm+');
        plot(lv_centroid(1), p_x(x_out), 'yo');
        plot(lv_centroid(1), p_x(x_in), 'go');


        subplot(no_of_rows, no_of_cols, 5);
        cla
        plot(p_x, p_y, 'b-');
        hold on;
        plot(p_x(x_last_trough), p_y(x_last_trough), 'ro');
        plot(p_x([1 end]), wall_thresh*[1 1],'c-');
        plot(p_x(x_last_peak), p_y(x_last_peak), 'rd');
        plot(p_x(x_out), p_y(x_out), 'gd');
        plot(p_x(x_in), p_y(x_in), 'gs');
        
        subplot(no_of_rows, no_of_cols, 6);
        hold on
        cm = paruly(no_of_frames);
        plot(ang, wall_thickness(i,:),'-','Color',cm(i,:));
        xlabel('Rotation Angle');
        ylabel('Wall thickness');

        drawnow;
%         pause(0.2);
        fr = getframe(gcf);
        writeVideo(v,fr);
        

    end
end

close(v)