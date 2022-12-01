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

%% Hist
imhist(cropped_Im)

%% binarization
Im_BW = imbinarize(cropped_Im,0.8);

figure
imshowpair(cropped_Im,Im_BW,'montage')

%% morphological operation
close all
se = strel('disk',2);
Dil = imdilate(Im_BW,se);

figure
imshowpair(cropped_Im,Dil,'montage')
title('Dilation')

se = strel('disk',9);
Op=imopen(Dil,se);

figure
imshowpair(cropped_Im,Op,'montage')
title('Opening')

%% Segmentation
%GT = uint8(GT);
GT = im2double(Op);
cropped_Im=im2double(cropped_Im);

%Seg = GT*cropped_Im;
Seg =  immultiply(cropped_Im, GT);
figure
imshow(Seg)
