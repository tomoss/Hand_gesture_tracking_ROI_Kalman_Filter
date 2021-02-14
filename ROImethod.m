clc; 
close all;
clear;

%% Init video

video_name = 'video2.mp4'; % name video: video1, video2, video3
vid = VideoReader(video_name); 
frameRate = vid.FrameRate;
nframes = vid.NumFrames; 
Height = vid.Height;
Width = vid.Width;

%% Init parameters - mean and standard deviation

if (strcmp(video_name,'video3.mp4') == 1)
    % parameters for video3
    u_cb = 120.3846; 
    u_cr = 150.7692;
    sigma_cb = 37.136041;
    sigma_cr = 13.80914;
else
    % parameters for video1, video2
    u_cb = 119.3846; 
    u_cr = 141.7692;
    sigma_cb = 8.136041;
    sigma_cr = 13.80914;
end

%% Display St image

RGB = read(vid,1); 
YCBCR = rgb2ycbcr(RGB);

Y  = YCBCR(:,:,1);
Cb = YCBCR(:,:,2);
Cr = YCBCR(:,:,3);

St = zeros(size(Cr,1), size(Cr,2));
for i = 1:size(Cr,1)
    for j = 1:size(Cr,2)
        if u_cr - sigma_cr < Cr(i,j) && Cr(i,j) < u_cr + sigma_cr 
            if u_cb - sigma_cb < Cb(i,j) && Cb(i,j) < u_cb + sigma_cb
                St(i,j) = 1;
            end
        end
    end
end

subplot(1,2,1); imshow(RGB); title('Imagine RGB');
subplot(1,2,2); imshow(St); title('Imagine St');

%% Display SM image

Ft1 = read(vid,1); Ft2 = read(vid,2);
Fdt = abs(Ft2 - Ft1);
R = Fdt(:, :, 1); G = Fdt(:, :, 2); B = Fdt(:, :, 3);
Fdg = 0.299*R + 0.587*G + 0.114*B; % Convert to gray

X = mean(Fdg, 'all');
T = 0.05 * X;

% Fdb = zeros(size(Fdg,1), size(Fdg,2));
% for i = 1:size(Fdg,1)
%     for j = 1:size(Fdg,2)
%         if double(Fdg(i,j))/255 >= T
%             Fdb(i,j) = 1;
%         end
%     end
% end
Fdb = imbinarize(Fdg,T);

figure,
subplot(2,2,1); imshow(Ft1); title('Image at t-1');
subplot(2,2,2); imshow(Ft2); title('Image at t');
subplot(2,2,3); imshow(Fdg); title('Fdg');
subplot(2,2,4); imshow(Fdb); title('Fdb');

SM = St & Fdb;
SM = medfilt2(SM);

figure,
subplot(1,2,1); imshow(RGB); title('RGB');
subplot(1,2,2); imshow(SM); title('SM');

%% Tracking using ROI based method

init = 0; 
for k=1:nframes-1
    
    im1 = read(vid,k);
    im2 = read(vid,k+1); 
    
    [FSM] = getHandPos(im1, im2, u_cb, u_cr, sigma_cb, sigma_cr);
    
    if (init == 0)
        [min_col, min_row, max_col, max_row, widthbbox, highbbox] = calculateROI(FSM);
        init = 1;
    end
    
    im1 = insertShape(im1,'Rectangle',[min_col min_row widthbbox highbbox],'LineWidth',2); % insereaza patratul pozitiei curente
    
    scenario_ROI(:,:,:,k) = im1;

    [min_col_left, min_col_right, min_row_top, right_widthbbox, left_widthbbox, top_highbbox] = getRegions(min_col, min_row, max_col, max_row, widthbbox, highbbox);
    
    right = imcrop(FSM,[min_col_right min_row right_widthbbox highbbox]);
    left  = imcrop(FSM,[min_col_left min_row left_widthbbox highbbox]);
    top   = imcrop(FSM,[min_col min_row_top widthbbox top_highbbox]);
    
    sumright = sum(right,'all');
    sumleft  = sum(left, 'all');
    sumtop   = sum(top,  'all');
    
    [min_right, max_left, min_top]= getMinMax(right, left, top);
    
    if ( sumright > sumleft ) % shift the ROI to the right
        right_displacement = right_widthbbox - min_right;
        new_min_col = min_col - right_displacement; 
        new_max_col = max_col - right_displacement;
        disp('Move right')
    elseif ( sumleft > sumright ) % shift the ROI to the left
        left_displacement = max_left;
        new_min_col = min_col + left_displacement; 
        new_max_col = max_col + left_displacement;
        disp('Move left')
    else
        new_min_col = min_col;
        new_max_col = max_col;
    end
    if ( sumtop > 0 )% shift the ROI either to top or bottom
        if ( min_top < 0.5 * top_highbbox) % shift to top
            new_min_row = min_row - ( 0.5 * top_highbbox - min_top );
            disp('Move top')
        elseif ( min_top > 0.5 * top_highbbox ) % shift to bottom
            new_min_row = min_row + ( min_top - 0.5 * top_highbbox );
            disp('Move bottom')
        end
    else
        disp('Not moving')
        new_min_row = min_row;
    end
    
    min_col = new_min_col;
    min_row = new_min_row;
    max_col = new_max_col;
    
end

%% Display and save video

implay(scenario_ROI, frameRate);

v = VideoWriter(strcat(video_name(1:6),'_ROI_output'), 'MPEG-4');
v.FrameRate = vid.FrameRate;
open(v);
writeVideo(v,scenario_ROI);
close(v);

