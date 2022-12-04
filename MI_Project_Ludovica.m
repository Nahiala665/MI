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
Im_d = Im(:,:,10);
Im_int = uint8(Im_d); % one slice 

figure
subplot(1, 2, 1), imshow(Im_d, [])
subplot(1, 2, 2), imshow(Im_int, []) %brighter

Im_bin = imbinarize(Im_d, 128);

figure
subplot(1,3,1), imshow(Im_bin,[]), title('Binarized')
subplot(1,3,2), imcontour(Im_int,1,'m')
subplot(1,3,3), imshow(Im_bin,[]), hold on, imcontour(Im_int,1,'m'), title('Contours')

%% cropping one slice 

[cropped_Im, d] = imcrop(Im_int);
close all

d = round(d); % i save the coords

figure
subplot(1,2,1), imshow(Im_int), title('Original')
subplot(1,2,2), imshow(cropped_Im), title('Cropped')

%% aumatic crop of all the slices 
% coordinate

v1 = d(1):d(1)+d(3);
v2 = d(2):d(2)+d(4);

Im_int_cropped = uint8(Im);
% figure
% montage(Im(v1, v2, :), cmap)

figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_int_cropped(v1,v2,i),[])
end

%% 
FT = fft2(Im_int_cropped);

IFT = ifft2(FT);


figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_int_cropped(v1,v2,i), [])
end
title('Original')

figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(IFT(v1,v2,i), [])
end
title('Reconstructed')

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
cropped_Im_d = im2double( Im_int_cropped );    

max_im_cropped = max(cropped_Im_d(:) );
min_im_cropped = min( cropped_Im_d(:) );

LOW_in = min_im_cropped;
HIGH_in = max_im_cropped;

LOW_out = 0;
HIGH_out = 1;
gamma = 5;
        % supposing linear transform
        % mapping input dynamic range into the whole dyamic range as output
        
for i = 1:20
    cropped_Im_modified(v1,v2,i) = imadjust(cropped_Im_d(v1,v2,i), [LOW_in HIGH_in], [LOW_out HIGH_out], gamma); % this works better
end

% Show the contrast change related to the adjustment operation

figure 
title('Im double modified gamma = 5')
for i = 1:20
    subplot(4, 5, i)
    imshow(cropped_Im_modified(v1,v2,i), [])
end


% negative image
%Im_neg = imcomplement( cropped_Im ); 

%figure
%imshow(Im_neg, [])
%title('NEGATIVE');

%cropped_Im = im2double( cropped_Im );  

%%%%%%%%%%%%%%%%%%%%
%% binarizing

Th = graythresh(cropped_Im_modified);      % Otsu's method

bin_cropped_Im = imbinarize(cropped_Im_modified, Th);

figure 
title('Im double modified gamma = 5')
for i = 1:20
    subplot(4, 5, i)
    imshow(bin_cropped_Im(v1,v2,i), [])
end

%%%%%%%%%%





%%%%%%%%%%
close all
for i = 1:20
    diff_cropped(v1, v2, i) = cropped_Im_d(v1, v2, i) - cropped_Im_modified(v1, v2, i);
end


figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(diff_cropped(v1,v2,i))
    colorbar
end
title('diff gamma = 4')


%% filtering
close all
% equalization of the histogram

h1 = histeq(cropped_Im, 256); % this works good
h2 = histeq(diff_cropped, 256);

% figure
% subplot(2, 2, 3)
% imhist(h1)
% subplot(2, 2, 4)
% imhist(h2)

figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(h1(v1,v2,i))
    colorbar
end
title('cropped Im histogram eq')

figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(h2(v1,v2,i))
    colorbar
end
title('cropped Im histogram diff eq')


%% good 
close all

h1_d = im2double( h1 );    

max_im_h1 = max( h1_d(:) );
min_im_h1 = min( h1_d(:) );

LOW_in = min_im_h1;
HIGH_in = max_im_h1;

LOW_out = 0;
HIGH_out = 1;
gamma = 5;
        % supposing linear transform
        % mapping input dynamic range into the whole dyamic range as output

for i = 1:20
    h1_modified(v1,v2,i) = imadjust(h1_d(v1,v2,i), [LOW_in HIGH_in], [LOW_out HIGH_out], gamma); % this works better
end

figure
for i = 1:20
    subplot(4, 5, i)
    imshow(h1_modified(v1,v2,i), [])
end
title('Im double modified gamma = 5')

% subplot(221), imshow(h1), title('h1 double')
% subplot(222), imshow(h1_modified), title('h1 double modified gamma = 5')
% subplot(223), imhist(h1), ylim([0 3000])
% subplot(224), imhist(h1_modified), ylim([0 3000])
%%
h1_bin = imbinarize(h1_modified); % the circle is enhanced

figure
for i = 1:20
    subplot(4, 5, i)
    imshow(h1_bin(v1,v2,i), [])
end

%% binarization
Im_BW = imbinarize(cropped_Im,0.8);

figure
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_BW(v1,v2,i), [])
end

%% morphological operation
close all

se = strel('disk',2);
Dil = imdilate(Im_BW,se);

figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Dil(v1,v2,i),[])
end
title('Dilation')

se = strel('disk',9);
Op=imopen(Dil,se);

figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Op(v1,v2,i),[])
end
title('Opening')


%% Segmentation
close all

%GT = uint8(GT);
GT = im2double(Op);
cropped_Im = im2double(cropped_Im);

%Seg = GT*cropped_Im;
Seg =  immultiply(cropped_Im, GT);

figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(Seg(v1,v2,i))
end
title('Segm')

%% video
close all 

figure
for i = 1:20
    imshow(Seg(v1,v2,i))
    hold on
    pause
end   



