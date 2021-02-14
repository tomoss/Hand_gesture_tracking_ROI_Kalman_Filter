function [SM]= getHandPos(im1, im2, u_cb, u_cr, sigma_cb, sigma_cr)
% Skin Color Segmentation

RGB = im1;
YCBCR = rgb2ycbcr(RGB);

Y =  YCBCR(:,:,1);
Cb = YCBCR(:,:,2);
Cr=  YCBCR(:,:,3);

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

%% Moving Object Segmentation

Ft1 = im1; Ft2 = im2;
Fdt = abs(Ft2 - Ft1);
R = Fdt(:, :, 1);
G = Fdt(:, :, 2);
B = Fdt(:, :, 3);
Fdg = 0.299*R + 0.587*G + 0.114*B; % Convert to gray

X = mean(Fdg, 'all');
T = 0.05 * X;

%Fdb = zeros(size(Fdg,1), size(Fdg,2));
% for i = 1:size(Fdg,1)
%     for j = 1:size(Fdg,2)
%         if double(Fdg(i,j))/255 >= T
%             Fdb(i,j) = 1;
%         end
%     end
% end
Fdb = imbinarize(Fdg,T); 

%% Defining The Hand Region

SM = St & Fdb;
SM = medfilt2(SM); 

end

