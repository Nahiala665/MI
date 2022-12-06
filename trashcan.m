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
