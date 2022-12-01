clear all
clc

%%
Image = load ('CMRIdata.mat'); % this is a struct

Im = Image.vol; % matrix of interest, this is a double

% double Gray scale
cmap = colormap('gray');
montage(Im, cmap) % to see all the slices --> 20 slices
title('double-GrayScale')
colorbar

Im_max = max(max(Im)); % to have the max value of intensity 

% we consider mid_slice 10
Im_d = Im(:,:,20);

Im_bin = imbinarize(Im_d, 128);

Im_int = uint8(Im_d);

figure
subplot(1,3,1), imshow(Im_bin,[]), title('Binarized')
subplot(1,3,2), imcontour(Im_int,1,'m')
subplot(1,3,3), imshow(Im_bin,[]), hold on, imcontour(Im_int,1,'m'), title('Contours')

%%
cropped_Im = imcrop(Im_int);

figure
subplot(1,2,1), imshow(Im_int), title('Original')
subplot(1,2,2), imshow(cropped_Im), title('Cropped')

% Th = graythresh(cropped_Im);      % Otsu's method
% 
% bin_cropped_Im = imbinarize(cropped_Im, Th);
% 
% figure
% subplot(1,2,1), imshow(cropped_Im, []), title('Original')
% subplot(1,2,2), imshow(bin_cropped_Im, []), title('logical-Binary')
% 
% FT = fft2(cropped_Im);
% 
% IFT = ifft2(FT);
% 
% figure
% subplot(1,2,1), imshow(cropped_Im, []), title('Original')
% subplot(1,2,2), imshow(IFT, []), title('Reconstructed')


%%
% negative image
Im_neg = imcomplement( cropped_Im ); 

figure
imshow(Im_neg, [])
title('NEGATIVE');

%cropped_Im = im2double( cropped_Im );  

%%

cropped_Im_d = im2double( cropped_Im );    

max_im_cropped = max( cropped_Im_d(:) );
min_im_cropped = min( cropped_Im_d(:) );

LOW_in = min_im_cropped;
HIGH_in = max_im_cropped;

LOW_out = 0;
HIGH_out = 1;
gamma = 4;
        % supposing linear transform
        % mapping input dynamic range into the whole dyamic range as output

cropped_Im_modified = imadjust(cropped_Im_d, [LOW_in HIGH_in], [LOW_out HIGH_out], gamma);

% Show the contrast change related to the adjustment operation
figure
subplot(221), imshow(cropped_Im_d)
subplot(222), imshow(cropped_Im_modified)
subplot(223), imhist(cropped_Im_d), ylim([0 3000])
subplot(224), imhist(cropped_Im_modified), ylim([0 3000])

diff_cropped = cropped_Im_d - cropped_Im_modified;

figure
subplot(211), imshow(diff_cropped), title('diff gamma = 4'), colorbar
subplot(212), imhist(diff_cropped), ylim([0 3000])

%% filtering

% h = fspecial('log',7,0.4);
% I2 = imfilter(cropped_Im_d,h);
% 
% figure
% imshow(I2)

%volumeViewer(Im)

% equalization of the histogram

h1 = histeq(cropped_Im, 256);
h2 = histeq(diff_cropped, 256);

figure
subplot(2, 2, 1)
imshow(h1)
subplot(2, 2, 2)
imshow(h2)
subplot(2, 2, 3)
imhist(h1)
subplot(2, 2, 4)
imhist(h2)

Im_bin1 = imbinarize(h2);

figure
imshow(Im_bin1,[]), title('Binarized')