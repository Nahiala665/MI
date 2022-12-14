clear all
clc

%%
close all
Images = load ('CMRIdata.mat'); % this is a struct

Im = Images.vol; % matrix of interest, this is a double

% double Gray scale
cmap = colormap('gray');
montage(Im, cmap) % to see all the slices --> 20 slices
title('double-GrayScale')
colorbar


%% We consider mid_slice 2
close all

Im_d = Im(:,:,2);
Im_int = uint8(Im_d); % one slice 

[cropped_Im, d] = imcrop(Im_int);

d = round(d); % we save the coords

figure
subplot(1,2,1), imshow(Im_int), title('Original')
subplot(1,2,2), imshow(cropped_Im), title('Cropped')

%% automatic crop of all the slices 
% coordinate

v1 = d(1):d(1)+d(3);
v2 = d(2):d(2)+d(4);

% Im_int_cropped = uint8(Im);

figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Im(v2,v1,i),[])
end

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


%% creating the disk for all (uncropped image size)
BG=size(Im_int);
LV = zeros(256,256,3,20);
figure
for i=1:20
    subplot(4,5,i)
    whiteImage = 0 * ones(BG(1), BG(2), 'uint8');
    center1 = center(i,1)+d(1);
    center2 = center(i,2)+d(2);
    J = insertShape(whiteImage,'filled-circle',[center1 center2 radius(i)],'color','white','opacity',1);
    LV(:,:,:,i)=im2double(J);
    imshow(J)
end

%% Segmentation (uncropped image)
figure
segmentation=zeros(256,256,20);
for i=1:20
    subplot(4,5,i)
    Full_im=im2double(Im_int_cropped(:,:,i));
    Crop_circle = LV(:,:,1,i);
    Seg = immultiply(Full_im, Crop_circle);
    segmentation(:,:,i)=Seg;
    imshow(Seg)
end

%% video of the segmented ROI
close all 

figure
for i = 1:20
    imshow(segmentation(:,:,i))
    hold on
    pause
end  

%% groundtruth

Im2 = Images.gsmask; % matrix of interest, this is a double

% double Gray scale
figure
for i = 1:20
    subplot(4,5,i)
    imshow(Im2(:,:,i));
end


%% Evaluate Image Segmentation Score
%Sensitivity (true positive rate) refers to the probability of a positive test, conditioned on truly being positive.
% Specificity (true negative rate) refers to the probability of a negative test, conditioned on truly being negative.
dice_index=zeros(1,20);
sensitivity_index = zeros(1,20);
specificity_index =zeros(1,20);

for i=1:20
    subplot(4,5,i)
    LV_BW = LV(:,:,:, i);
    LV_BW=imbinarize(LV_BW(:, :, 1));
    GT_BW=imbinarize(Im2(:,:,i));
    [sensitivity_index(1,i),specificity_index(1,i),dice_index(1,i)] = SegmentationPerformance(GT_BW,LV_BW);
    imshowpair(LV_BW, GT_BW)
    sensitivity=round(sensitivity_index(1,i),3);
    specificity = round(specificity_index(1,i),3);
    similarity = round(dice_index(1,i),3);
    title(['slice' num2str(i)])
    txt = {['Dice : ' num2str(similarity)], ['TPR : ' num2str(sensitivity)],['TNR : ' num2str(specificity)]};

    text(0,210,txt,'FontSize',8,'Color','white')
end

mean_dice=mean(dice_index);