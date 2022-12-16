clear all
clc

%% Pre-explanation

% Our segmentation method will require the intervention of the user at 2 moments of the code. 

% The first time, it will require the user to crop a rectangle around the left ventricle of the heart on one of the slice only
% The second time, it will require the user to draw the diameter of the inner circle of the left ventricle on one slice only

% The rest of the code is fully automated and use a circle recognition method

%% Loading of the images
close all
Images = load ('CMRIdata.mat'); 

Im = Images.vol; % matrix of interest

% double Gray scale
figure
cmap = colormap('gray');
montage(Im, cmap) % to see the 20 slices
title('double-GrayScale')
colorbar


%% We crop on one slice
close all

Im_int = uint8(Im(:,:,2)); % we select one slice 


% We will crop a rectangle around the limit of the left ventricule
[cropped_Im, d] = imcrop(Im_int);

d = round(d); % we save the coordinates of our cropping

figure
subplot(1,2,1), imshow(Im_int), title('Original')
subplot(1,2,2), imshow(cropped_Im), title('Cropped')

%% automatic crop of all the slices 
% coordinates
v1 = d(1):d(1)+d(3);
v2 = d(2):d(2)+d(4);


figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Im(v2,v1,i),[])
end


%% Diameter of the circle
% Here, the user is asked to select the diameter inner circle of the left ventricule 
% (the white part only, not the darker wall around)


close all
Im_double = im2double( Im ); 
imshow(Im_double(v2,v1,20),[])
diam = drawline;
pos = diam.Position;
diffPos = diff(pos);
diameter = hypot(diffPos(1),diffPos(2));


%% Automatic circle recognition
% With the diameter of the circle given, we use the function imfindcircle to recognize the shape 
% of the left ventricle on the images

low = round(diameter/2) - 4;
up = round(diameter/2) + 4;

center=zeros(20,2);
radius=zeros(20,1);

figure
for i = 1:20
    subplot(4,5,i)
    [centers,radii] = imfindcircles(Im_double(v2,v1,i),[low up],'Sensitivity',0.95); %explain in the paper
    center(i,:)=centers(1,:);
    radius(i)=radii(1);
    imshow(Im_double(v2,v1,i),[])
    h = viscircles(centers,radii,'LineStyle','--');
end


%% creating the disk 
% Here, we create our proposition of segmentation shape for the left ventricule. 
% For that, we draw circle of the same diameter and the same center than the circle recognized in the previous part.
% Later, we will compare this shape proposition with the groundTruth to see how efficient our method is

BG=size(Im_int); %size of the background
LV = zeros(256,256,3,20); % left ventricule shape as detected before
figure
for i=1:20
    subplot(4,5,i)
    whiteImage = 0 * ones(BG(1), BG(2), 'uint8'); %we create a black background
    center1 = center(i,1)+d(1); %we add d(1) and d(2) to move the center on our full sized image
    center2 = center(i,2)+d(2);
    J = insertShape(whiteImage,'filled-circle',[center1 center2 radius(i)],'color','white','opacity',1); %white disks
    LV(:,:,:,i)=im2double(J);
    imshow(J)
end

%% Segmentation 
% In this part, we segmentate our image with the disk shapes.

figure
segmentation=zeros(256,256,20);
for i=1:20
    subplot(4,5,i)
    Full_im=Im_double(:,:,i);
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
% We now load the groundtruth given in our dataset

Im2 = Images.gsmask; 

% double Gray scale
figure
for i = 1:20
    subplot(4,5,i)
    imshow(Im2(:,:,i));
end

%% True segmentation
% We segmentate our images with the groundtruth shape

segmentationTruth=zeros(256,256,20);
Im_double_Truth = im2double(Im2);

figure
for i=1:20
    subplot(4,5,i)
    Full_im=Im_double(:,:,i);
    Crop_circle_Truth = Im_double_Truth(:,:,i);
    SegTruth = immultiply(Full_im, Crop_circle_Truth);
    segmentationTruth(:,:,i)=SegTruth;
    imshow(SegTruth)
end


%% Evaluate Image Segmentation Score
% Here, we will compare our segmentation shape with the groundtruth. 
% To do that, we created a function called SegmentationPerformance that
% calculates the sensitivity, the specificity and the dice index of our methods, comparing the groundtruth and the segmentation shape 

% Sensitivity (true positive rate) refers to the probability of a positive test, conditioned on truly being positive.
% Specificity (true negative rate) refers to the probability of a negative test, conditioned on truly being negative.
% Dice index measures the similarity between two sets of data

% empty matrix to store our results
dice_index=zeros(1,20);
sensitivity_index = zeros(1,20);
specificity_index =zeros(1,20);

% the segmentation shape and grountruth must be binary images
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

% mean of our segmentation scores
mean_dice=round(mean(dice_index),3);
mean_sensitivity = round(mean(sensitivity_index),3);
mean_specificity=round(mean(specificity_index),3);

%% Plotting the segmentation score over the time frames

figure
plot(dice_index(1,:))
ylim([0 1])
grid on
hold on
plot(sensitivity_index(1,:),'Color','red')
plot(specificity_index(1,:),'Color','blue')
legend('similarity', 'sensitivity','specificity')
title('Segmentation evaluation over time frames')
    txt = {['Dice mean: ' num2str(mean_dice)], ['TPR mean : ' num2str(mean_sensitivity)],['TNR mean : ' num2str(mean_specificity)]};
    text(2,0.2,txt,'FontSize',10,'Color','black')


%% Quantification of the cross-sectional volume of the left ventricle

Volume = zeros(20,2);

% With our segmentation methods
for i=1:20
    LV_BW = LV(:,:,:, i);
    LV_BW=imbinarize(LV_BW(:, :, 1));
    Volume(i,1)=sum(LV_BW(:));
end

% GroundTruth Volume
for i=1:20
    GT_BW=imbinarize(Im2(:,:,i));
    Volume(i,2)=sum(GT_BW(:));
end

% Ploting the volumes obtained
figure
plot(Volume(:,1))
hold on
plot(Volume(:,2),'Color','red')
legend('calculated volume', 'real volume')
title('Comparison of the calculated and true volume of the left ventricule over time frames')

%% video of the segmented ROI with volume
close all 

figure
for i = 1:20
    imshow(segmentation(:,:,i))
    vol1 = Volume(i,1);
    vol2 = Volume(i,2);
    title(['Calculated volume :' num2str(vol1)])
    txt = {['True volume : ' num2str(vol2)]};
    text(100,0,txt,'FontSize',10,'Color','white')
    hold on
    pause
end  



