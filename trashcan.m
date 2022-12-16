%% Intensity
Im_max = max(max(Im)); % to have the max value of intensity -> do you use that later ? 

%% contour
figure
subplot(1,3,1), imshow(Im_bin,[]), title('Binarized')
subplot(1,3,2), imcontour(Im_int,1,'m')
subplot(1,3,3), imshow(Im_bin,[]), hold on, imcontour(Im_int,1,'m'), title('Contours')
% do we use the contours ? 

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

%% binarizing
Th = graythresh(diff_cropped); 
Im_BW = imbinarize(diff_cropped,Th);

figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_BW(v2,v1,i), [])
end

%% morphological operation
close all

se = strel('disk',2);
Clo = imclose(Im_BW,se);

figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Clo(v2,v1,i),[])
end
title('Closing')

% se = strel('disk',1);
% Ero = imerode(Im_BW,se);
% 
% figure()
% for i = 1:20
%     subplot(4, 5, i)
%     imshow(Ero(v2,v1,i),[])
% end
% title('Dilation')

% se = strel('disk',4);
% Op=imclose(Dil,se);
% 
% figure()
% for i = 1:20
%     subplot(4, 5, i)
%     imshow(Op(v2,v1,i),[])
% end
% title('Opening')

%% Create the disk
disk10 = circle(center(10,1),center(10,2),radius(10));

disk10=uint8(disk10);

figure
imshow(disk10)

% Seg =  immultiply(Im_int, disk10);
% figure
% imshow(Seg);

%% legal trick
% J = uint8(cropped_Im);
J=imadjust(cropped_Im,[0 0.01]);

figure
imshow(cropped_Im)
figure
imshow(J)
I = insertShape(J,'filled-circle',[center(10,1) center(10,2) radius(10)],'color',[1 1 1],'opacity',1); 
figure
imshow(I)

%% automatic legal trick

figure
for i=1:20
    subplot(4,5,i)
    J=imadjust(Im_int_cropped(v2,v1,i),[0 0.0001]);
    I = insertShape(J,'filled-circle',[center(i,1) center(i,2) radius(i)],'color',[1 1 1],'opacity',1);
    Seg = immultiply(Im_int_cropped(v2,v1,i),I(:,:,1));
    imshow(Seg)
end

%% save image
Im_d = Im(:,:,20);
Im_int = uint8(Im_d);
imwrite(Im_int, "im20.png")

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


%% dice similarity
dice_index=zeros(1,20);
figure
for i=1:20
    subplot(4,5,i)
    LV_BW = LV(:,:,:, i);
    LV_BW=imbinarize(LV_BW(:, :, 1));
    GT_BW=imbinarize(GT(:,:,1,i));
    similarity = dice(GT_BW,LV_BW);
    dice_index(1,i)=similarity;
    imshowpair(LV_BW, GT_BW)
    title(['dice index : ' num2str(similarity)])
end

mean_dice=mean(dice_index);


%% groundtruth (done by illustrator)
GT=zeros(256,256,3,5);

% extract the groundtruth
for i=1:20
    jpgFilename = sprintf('GT%d.png', i);
    GT(:,:,:,i)=imread(jpgFilename);
end

%% We consider mid_slice 10
Im_d = Im(:,:,2);
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

%% Circle recognition for slide 10
close all

low = round(diameter/2) - 5;
up = round(diameter/2) + 5;
[centers,radii] = imfindcircles(Image(v2,v1,20),[low up],'Sensitivity',0.97);

imshow(Image(v2,v1,20))
h = viscircles(centers,radii,'LineStyle','--');


%% creating the disk for sl0
S=size(cropped_Im);
whiteImage = 255 * ones(S(1), S(2), 'uint8');
I = insertShape(whiteImage,'filled-circle',[center(10,1) center(10,2) radius(10)],'color',[1 1 1],'opacity',1); 
figure
imshow(I)

%% Segmentation (on the cropped image
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
