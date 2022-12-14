clear all
clc

%%
close all
Image = load ('CMRIdata.mat'); % this is a struct

Im = Image.vol; % matrix of interest, this is a double

% double Gray scale
cmap = colormap('gray');
montage(Im, cmap) % to see all the slices --> 20 slices
title('double-GrayScale')
colorbar


%% We consider mid_slice 10
Im_d = Im(:,:,10);
Im_int = uint8(Im_d); % one slice 
% we need to convert to uint8 to use the imcrop

figure
subplot(1, 2, 1), imshow(Im_d, [])
subplot(1, 2, 2), imshow(Im_int, []) %brighter

Im_bin = imbinarize(Im_d, 128);

%% cropping one slice 
% we should crop from the slice where the LV is the biggest
[cropped_Im, d] = imcrop(Im_int);
close all

d = round(d); % we save the coords

figure
subplot(1,2,1), imshow(Im_int), title('Original')
subplot(1,2,2), imshow(cropped_Im), title('Cropped')

%% automatic crop of all the slices 
% coordinate

v1 = d(1):d(1)+d(3);
v2 = d(2):d(2)+d(4);

Im_int_cropped = uint8(Im);

figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_int_cropped(v2,v1,i),[])
end

%% Resolution work
close all

FT = fft2(Im_int_cropped);

IFT = ifft2(FT);


figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_int_cropped(v2,v1,i), [])
end
title('Original')

figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(IFT(v2,v1,i), [])
end
title('Reconstructed')

%% make it darker
cropped_Im_d = im2double( Im_int_cropped );    

% max_im_cropped = max(cropped_Im_d(:) );
% min_im_cropped = min( cropped_Im_d(:) );
% 
% LOW_in = min_im_cropped;
% HIGH_in = max_im_cropped;
% 
% LOW_out = 0;
% HIGH_out = 1;

gamma = 5;
        % supposing linear transform
        % mapping input dynamic range into the whole dyamic range as output
        
for i = 1:20
    cropped_Im_modified(v2,v1,i) = imadjust(cropped_Im_d(v2,v1,i), [0 1], [0 1], gamma); % this works better
end

% Show the contrast change related to the adjustment operation

figure 
title('Im double modified gamma = 5')
for i = 1:20
    subplot(4, 5, i)
    imshow(cropped_Im_modified(v2,v1,i), [])
end

%% difference
close all 

for i = 1:20
    diff_cropped(v2, v1, i) = cropped_Im_d(v2, v1, i) - cropped_Im_modified(v2, v1, i);
end


figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(diff_cropped(v2,v1,i))
    colorbar
end
title('diff gamma = 4')

%% Diameter of the circle
close all
Image = im2double(diff_cropped);
imshow(Image(v2,v1,20))
diam = drawline;
pos = diam.Position;
diffPos = diff(pos);
diameter = hypot(diffPos(1),diffPos(2));

%% Circle recognition for slide 10
close all

low = round(diameter/2) - 5;
up = round(diameter/2) + 5;
[centers,radii] = imfindcircles(Image(v2,v1,20),[low up],'Sensitivity',0.97);

imshow(Image(v2,v1,20))
h = viscircles(centers,radii,'LineStyle','--');

%% Automatic circle recognition
figure
center=zeros(20,2);
radius=zeros(20,1);
for i = 1:20
    subplot(4,5,i)
    [centers,radii] = imfindcircles(Image(v2,v1,i),[low up],'Sensitivity',0.97);
    center(i,:)=centers(1,:);
    radius(i)=radii(1);
    imshow(Image(v2,v1,i))
    h = viscircles(centers,radii,'LineStyle','--');
end

%% creating the disk for sl0
S=size(cropped_Im);
whiteImage = 255 * ones(S(1), S(2), 'uint8');
I = insertShape(whiteImage,'filled-circle',[center(10,1) center(10,2) radius(10)],'color',[1 1 1],'opacity',1); 
figure
imshow(I)

%% cropping
Seg =  immultiply(cropped_Im, I(:,:,1));
figure
imshow(Seg)

%% automatic legal trick
S=size(Im_int_cropped(v2,v1,1));
figure
for i=1:20
    subplot(4,5,i)
    whiteImage = 255 * ones(S(1), S(2), 'uint8');
    I = insertShape(whiteImage,'filled-circle',[center(i,1) center(i,2) radius(i)],'color',[1 1 1],'opacity',1);
    Seg = immultiply(Im_int_cropped(v2,v1,i),I(:,:,1));
    imshow(Seg)
end

%% creating the disk for all
BG=size(Im_int);
LV = zeros(256,256,3,20);
figure
for i=1:20
    subplot(4,5,i)
    whiteImage = 255 * ones(BG(1), BG(2), 'uint8');
    center1 = center(i,1)+d(1);
    center2 = center(i,2)+d(2);
    J = insertShape(whiteImage,'filled-circle',[center1 center2 radius(i)],'color',[1 1 1],'opacity',1);
    LV(:,:,:,i)=im2double(J);
    imshow(J)
end

%% import the groundtruh for the slice 1
GT1 = imread ('GT1.png');
GT1=imbinarize(GT1(:,:,1));

figure
imshow(GT1)

LV1 = LV(:,:,:, 1);
LV1=imbinarize(LV1(:, :, 1));

figure
imshow(LV1);

similarity = dice(LV1,GT1);
figure
imshowpair(LV1, GT1)
title(['dice index : ' num2str(similarity)])

%% for 5 slices
GT=zeros(256,256,3,5);

% extract the groundtruth
for i=1:5
    jpgFilename = sprintf('GT%d.png', i);
    GT(:,:,:,i)=imread(jpgFilename);
end

% dice similarity
figure
for i=1:5
    subplot(1,5,i)
    LV_BW = LV(:,:,:, i);
    LV_BW=imbinarize(LV_BW(:, :, 1));
    GT_BW=imbinarize(GT(:,:,1,i));
    similarity = dice(LV_BW,GT_BW);
    imshowpair(LV_BW, GT_BW)
    title(['dice index : ' num2str(similarity)])
end

