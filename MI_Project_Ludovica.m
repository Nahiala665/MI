clear all
clc

%%
Image = load ('CMRIdata.mat'); % this is a struct

Im = Image.vol; % matrix of interest, this is a double

% double Gray scale
figure
cmap = colormap('gray');
montage(Im, cmap) % to see all the slices --> 20 slices
title('double-GrayScale')
%colorbar

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
close all

[cropped_Im, d] = imcrop(Im_int);
close all

d = round(d); % we save the coords

figure
subplot(1,2,1), imshow(Im_int), title('Original')
subplot(1,2,2), imshow(cropped_Im), title('Cropped')

%% automatic crop of all the slices 
% coordinate
close all

v1 = d(1):d(1)+d(3);
v2 = d(2):d(2)+d(4);

Im_int_cropped = uint8(Im);

figure()
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_int_cropped(v2,v1,i),[])
end
sgtitle('All slices cropped')
figure
for i = 1:20
    subplot(4, 5, i)
    imshow(Im_int_cropped(:,:,i),[])
end


%% make it darker
close all

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

for i = 1:20
    subplot(4, 5, i)
    imshow(cropped_Im_modified(v2,v1,i), [])
end
sgtitle('Im double modified gamma = 5')
%% difference
close all 

for i = 1:20
    diff_cropped(v2, v1, i) = cropped_Im_d(v2, v1, i) - cropped_Im_modified(v2, v1, i);
end


figure 
for i = 1:20
    subplot(4, 5, i)
    imshow(diff_cropped(v2,v1,i))
    %colorbar
end
sgtitle('Diff gamma = 5')

%% Diameter of the circle
close all
Image = im2double(diff_cropped);
imshow(Image(v2,v1,20))
d = drawline;
pos = d.Position;
diffPos = diff(pos);
diameter = hypot(diffPos(1),diffPos(2));

%% Circle recognition
close all

low = round(diameter/2) - 5;
up = round(diameter/2) + 5;
[centers,radii] = imfindcircles(Image(v2,v1,20),[low up],'Sensitivity',0.97);

imshow(Image(v2,v1,20))
h = viscircles(centers,radii,'LineStyle','--');

%% Automatic circle recognition
close all

figure
center=zeros(20,2);
radius=zeros(20,1);
for i = 1:20
    subplot(4,5,i)
    [centers,radii] = imfindcircles(Image(v2,v1,i),[low up],'Sensitivity',0.96);
    center(i,:)=centers(1,:);
    radius(i)=radii(1);
    imshow(Image(v2,v1,i))
    h = viscircles(centers,radii,'LineStyle','--');
end



%% legal trick
% J = uint8(cropped_Im);
close all

J=imadjust(cropped_Im,[0 0.01]);

figure
subplot(1, 3, 1)
imshow(cropped_Im), title('Cropped Image')
subplot(1, 3, 2)
imshow(J), title('White Image for the background')

I = insertShape(J,'FilledCircle',[center(10,1) center(10,2) radius(10)],'color',[1 1 1],'opacity',1); 

subplot(1, 3, 3)
imshow(I), title('Black circle in the ROI')


%% cropping
Seg =  immultiply(cropped_Im, I(:,:,1));
figure
imshow(Seg), title('Segmentation of the ROI, white backgroung')


%% creating the disk for sl0
S=size(cropped_Im);
whiteImage = 255 * ones(S(1), S(2), 'uint8');
I = insertShape(whiteImage,'FilledCircle',[center(10,1) center(10,2) radius(10)],'color',[1 1 1],'opacity',1); 
figure
imshow(I)

%% cropping
Seg =  immultiply(cropped_Im, I(:,:,1));
figure
imshow(Seg)
title('Segmentation of the midslice')

%% automatic legal trick
close all 

S=size(Im_int_cropped(v2,v1,1));
figure
for i=1:20
    subplot(4,5,i)
    whiteImage = 0 * ones(S(1), S(2), 'uint8'); 
    I = insertShape(whiteImage,'FilledCircle',[center(i,1) center(i,2) radius(i)],'color',[1 1 1],'opacity',1);
    Seg = immultiply(Im_int_cropped(v2,v1,i),I(:,:,1));
    imshow(Seg)
end

sgtitle('Segmentation of all the slices')

%% video of the segmented ROI
close all 
   
figure
for i = 1:20
    Seg = immultiply(Im_int_cropped(v2,v1,i),I(:,:,1));
    imshow(Seg)
    hold on
    pause
end  


%%
%%%%%%%


%% Circle recognition
close all

low = round(diameter/2) - 5;
up = round(diameter/2) + 5;
[centers,radii] = imfindcircles(Im_int_cropped(:,:,20),[low up],'Sensitivity',0.96);

imshow(Im_int_cropped(:,:,20))
h = viscircles(centers,radii,'LineStyle','--');

%% Automatic circle recognition
close all

figure
center=zeros(20,2);
radius=zeros(20,1);
for i = 1:20
    subplot(4,5,i)
    [centers,radii] = imfindcircles(Im_int_cropped(:,:,20),[low up],'Sensitivity',0.96);
    center(i,:)=centers(2,:);
    radius(i)=radii(2);
    imshow(Im_int_cropped(:,:,20))
    h = viscircles(centers,radii,'LineStyle','--');
end

%%
S1 = size(Im_int_cropped); % size of the entire image 
whiteImage = 255 * ones(S1(1), S1(2), 'uint8'); 
I = insertShape(whiteImage,'FilledCircle',[center(10,1) center(10,2) radius(10)],'color',[1 1 1],'opacity',1); 
% figure
% imshow(whiteImage, [])

figure
imshow(I)

%% cropping
Seg =  immultiply(Im_int_cropped(:,:,1), I(:,:,1));
figure
imshow(Seg)
title('Segmentation of the midslice')

%% automatic legal trick
close all 

S=size(Im_int_cropped(:,:,1));
figure
for i=1:20
    subplot(4,5,i)
    whiteImage = 0 * ones(S(1), S(2), 'uint8'); 
    I = insertShape(whiteImage,'FilledCircle',[center(i,1) center(i,2) radius(i)],'color',[1 1 1],'opacity',1);
    Seg = immultiply(Im_int_cropped(:,:,i),I(:,:,1));
    imshow(Seg)
end

sgtitle('Segmentation of all the slices')

%% video of the segmented ROI
close all 

figure
for i = 1:20
    Seg = immultiply(Im_int_cropped(:,:,i),I(:,:,1));
    imshow(Seg)
    hold on
    pause
end  
