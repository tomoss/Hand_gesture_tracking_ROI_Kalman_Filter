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

%% Parameters for Kalman Filter
dt=1;
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0 ; 0 0 0 1];
H = [1 0 0 0; 0 1 0 0];
B=0;
u = 0;
vt = 0;
Q=eye(4);
R = eye(2);
Pt_1 = eye(size(A,1));
f = 0.1845;
Ta = 1; % acceleration threshold

%% Tracking using Kalman

% Init
im1 = read(vid,1); im2 = read(vid,2);
FSM = getHandPos(im1, im2, u_cb, u_cr, sigma_cb, sigma_cr);
[min_col_1, min_row_1, max_col, max_row, widthbbox, highbbox] = calculateROI(FSM); 

im1 = read(vid,2); im2 = read(vid,3);
FSM = getHandPos(im1, im2, u_cb, u_cr, sigma_cb, sigma_cr);
[min_col, min_row, ~, ~, ~, ~] = calculateROI(FSM); 

xt_1 = [min_col; min_row; (min_col - min_col_1)/dt; (min_row - min_row_1)/dt];
xtp_1 = xt_1;

for k=3:nframes-1
    
    im1 = read(vid,k); im2 = read(vid,k+1);
    FSM = getHandPos(im1, im2, u_cb, u_cr, sigma_cb, sigma_cr);
    
    Ptp = A*Pt_1*A' + Q;
    Kt = Ptp*H'*inv(H*Ptp*H'+R);

    xt = A*xt_1 + B*u;
    zt = H*xt + vt;

    disp('Before correction')
    xtp = A*xtp_1 + B*u
     
    disp('After correction')
    xtp = xtp + Kt*(zt - H*xtp)

    % Estimate state
    im1 = insertShape(im1,'Rectangle',[xtp(1) xtp(2) widthbbox highbbox],'LineWidth',2,'Color','red');
    scenario_Kalman(:,:,:,k) = im1;

    Pt = (eye(size(Pt_1,1)) - Kt*H)*Ptp;

    % Adjust Kalman
    ax = (xtp(3) - xtp_1(3))/dt;
    ay = (xtp(4) - xtp_1(4))/dt;
    if( abs(ax) >= Ta || abs(ay) >= Ta)
        Wq=0.75/Ta;
        Wr=7.5/Ta/Ta;
    else
        Wq=0.25/Ta;
        Wr=7.5/Ta;
    end
    Q=Wq * eye(4);
    R=Wr*[f f/40; f/40 f];

    % Update Scanning Regions
    [min_col_left, min_col_right, min_row_top, right_widthbbox, left_widthbbox, top_highbbox] = getRegions(xtp(1), xtp(2), max_col, max_row, widthbbox, highbbox);

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
    elseif ( sumleft > sumright ) % shift the ROI to the left
        left_displacement = max_left;
        new_min_col = min_col + left_displacement; 
        new_max_col = max_col + left_displacement;
    else
        new_min_col = min_col;
        new_max_col = max_col;
    end
    if ( sumtop > 0 )% shift the ROI either to top or bottom
        if ( min_top < 0.5 * top_highbbox) % shift to top
            new_min_row = min_row - ( 0.5 * top_highbbox - min_top );
        elseif ( min_top > 0.5 * top_highbbox ) % shift to bottom
            new_min_row = min_row + ( min_top - 0.5 * top_highbbox );
        end
    else
        new_min_row = min_row;
    end

    min_col = new_min_col;
    min_row = new_min_row;
    max_col = new_max_col;

    vt = [min_col - xtp(1); min_row - xtp(2)];

    Pt_1 = Pt;
    xtp_1 = xtp;
    xt_1 = xtp;
 
end

%% Display and save video

implay(scenario_Kalman, frameRate);

v = VideoWriter(strcat(video_name(1:6),'_Kalman_output'), 'MPEG-4');
v.FrameRate = vid.FrameRate;
open(v);
writeVideo(v,scenario_Kalman);
close(v);








